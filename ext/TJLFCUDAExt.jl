module TJLFCUDAExt

using CUDA
using TJLF
using LinearAlgebra
import Random

# Per-device concurrency gate. A single ReentrantLock per GPU serialized every CUSOLVER dispatch.
# The A100 sits idle between one-at-a-time small solves, which is why MPS teams (separate
# processes = separate CUDA contexts) were needed for 2-4x scaling.
#
# EMPIRICAL FINDING (build/ad/batch_gpu_slots_validate_gpu.sh + smoke test): running the eigenvalue
# path (getrf/getrs + Xgeev) concurrently from multiple Julia tasks on ONE shared context corrupts
# results — a bitwise smoke test with 8 slots returned wrong eigenvalues on ~2% of 300 concurrent
# solves (max abs diff O(10)). CUDA.jl's per-task handles are NOT sufficient; Xgeev shares
# non-reentrant device state within a context. Therefore in-process concurrency is a correctness
# dead-end for the eigensolve. The only safe ways to fill the GPU are: separate contexts (MPS), or
# a single BATCHED cuSOLVER call over many matrices (the batched wave API — see TJLF._gpu_*_batch!).
#
# `TJLF_GPU_SLOTS` (default 1) sets concurrent solves per device:
#   * 1  -> Semaphore(1) == exact mutual exclusion == legacy per-device lock. THE ONLY SAFE VALUE
#           for the eigenvalue path; bit-for-bit identical dispatch behavior.
#   * >1 -> overlap N solves per device. UNSAFE for eigenvalues (see above); retained only for
#           experiments on the LU-only path. A one-time warning is emitted.
const _SLOTS_GUARD = ReentrantLock()
const _GPU_DEVICE_SEMS = Dict{Int,Base.Semaphore}()
const _GPU_DEVICE_WARMED = Dict{Int,Bool}()
const _SLOTS_WARNED = Threads.Atomic{Bool}(false)

function _gpu_slots()
    n = max(1, parse(Int, get(ENV, "TJLF_GPU_SLOTS", "1")))
    if n > 1 && !Threads.atomic_xchg!(_SLOTS_WARNED, true)
        @warn "TJLF_GPU_SLOTS=$n > 1: concurrent per-context eigensolves are NOT bitwise-safe " *
              "(cuSOLVER Xgeev corrupts results). Use MPS teams or the batched API for real " *
              "concurrency; slots>1 is for experiments on the LU-only path only."
    end
    return n
end

function _device_sem(id::Integer)
    lock(_SLOTS_GUARD) do
        get!(() -> Base.Semaphore(_gpu_slots()), _GPU_DEVICE_SEMS, Int(id))
    end
end

# Serialize the FIRST CUSOLVER touch on a device so concurrent tasks (slots>1) can't race
# cusolverDnCreate. `warm` runs a trivial getrf under the exclusive guard exactly once per device;
# subsequent calls short-circuit lock-free on the populated flag.
function _warmup_device!(id::Integer)
    get(_GPU_DEVICE_WARMED, Int(id), false) && return
    lock(_SLOTS_GUARD) do
        get(_GPU_DEVICE_WARMED, Int(id), false) && return
        z = CUDA.zeros(ComplexF64, 1, 1)
        ip = CUDA.CuArray{Int32}(undef, 1)
        try
            CUDA.CUSOLVER.getrf!(z, ip)
        catch
            # a failed 1x1 warmup must not wedge the device; real solves will surface errors
        end
        _GPU_DEVICE_WARMED[Int(id)] = true
    end
end

# Pick a device (round-robin across visible GPUs), warm it, then run `f` while holding one of the
# device's concurrency slots. With slots=1 this is mutual exclusion identical to the old lock.
function _with_device_slot(f)
    devs = _devices()
    dev = devs[mod1(Threads.atomic_add!(_GPU_RR, 1) + 1, length(devs))]
    CUDA.device!(dev)
    id = CUDA.deviceid(dev)
    _warmup_device!(id)
    sem = _device_sem(id)
    Base.acquire(sem)
    try
        return f()
    finally
        Base.release(sem)
    end
end

# Round-robin counter + cached device list to spread successive solves across all
# visible GPUs. Lazily initialised (double-checked) so it works whether or not
# CUDA is functional at load time.
const _GPU_RR = Threads.Atomic{Int}(0)
const _DEVICES = Ref{Vector{CuDevice}}(CuDevice[])
function _devices()
    devs = _DEVICES[]
    isempty(devs) || return devs
    lock(_SLOTS_GUARD) do
        devs = _DEVICES[]
        isempty(devs) || return devs
        devs = collect(CUDA.devices())
        _DEVICES[] = devs
        return devs
    end
end

# ── Batched shift-invert dominant-mode eigensolver (see TJLF/src/tjlf_batched_si.jl) ──────────
# Every step is a batched cuBLAS primitive over all P pencils at once, so it is concurrency-safe
# (unlike Xgeev) and keeps the A100 busy. Per shift σ: G = A-σB, batched LU (getrf), then the
# subspace operator S x = U⁻¹L⁻¹P(Bx) applied to the M-column block via batched GEMM + pivot
# gather + batched triangular solves; Q steps of Cholesky-QR subspace iteration keep only the
# M×M Gram matrices on the host; a tiny host M×M eig gives the Ritz μ → λ = σ + 1/μ.

# dst[i,c,p] = src[perm[i,p],c,p]  (apply per-pencil LU row permutation to the M-block)
function _si_gather_kernel!(dst, src, perm, n)
    i = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    c = CUDA.blockIdx().y
    p = CUDA.blockIdx().z
    if i <= n
        @inbounds dst[i, c, p] = src[perm[i, p], c, p]
    end
    return
end
function _si_gather_rows!(dst, src, perm)
    n, mm, P = size(dst)
    CUDA.@cuda threads=256 blocks=(cld(n, 256), mm, P) _si_gather_kernel!(dst, src, perm, n)
    return dst
end

# LAPACK ipiv (sequential swaps) -> permutation with (P·b)[i] = b[perm[i]]. Host, O(n) per pencil.
function _si_ipiv_to_perm(ipiv::AbstractMatrix{<:Integer})
    n, P = size(ipiv)
    perm = Matrix{Int32}(undef, n, P)
    idx = Vector{Int}(undef, n)
    for p in 1:P
        @inbounds for i in 1:n; idx[i] = i; end
        @inbounds for i in 1:n
            j = ipiv[i, p]
            j != i && ((idx[i], idx[j]) = (idx[j], idx[i]))
        end
        @inbounds for i in 1:n; perm[i, p] = idx[i]; end
    end
    return perm
end

# GPU Cholesky-QR: X <- X R⁻¹ with R = chol(XᴴX). Only the M×M Gram/Rinv cross PCIe.
function _si_cholqr!(X, gram, rinv, tmp)
    _, mm, P = size(X)
    CUDA.CUBLAS.gemm_strided_batched!('C','N', ComplexF64(1), X, X, ComplexF64(0), gram)
    Gh = Array(gram); Rinv = similar(Gh)
    Threads.@threads for p in 1:P
        @views begin
            H = Hermitian((Gh[:, :, p] .+ Gh[:, :, p]') ./ 2)
            Rinv[:, :, p] .= inv(cholesky(H).U)
        end
    end
    copyto!(rinv, Rinv)
    CUDA.CUBLAS.gemm_strided_batched!('N','N', ComplexF64(1), X, rinv, ComplexF64(0), tmp)
    copyto!(X, tmp)
    return X
end

function _si_shift!(cands, A3, B3, σ, X0h, bufs, M, Q, mu_tol)
    n, _, P = size(A3)
    X, Y, permbuf, tmp, gram, rinv = bufs
    copyto!(X, repeat(reshape(X0h, n, M, 1), 1, 1, P))
    G = A3 .- σ .* B3
    piv, _, _ = CUDA.CUBLAS.getrf_strided_batched!(G, true)   # LU in place
    perm = CuArray(_si_ipiv_to_perm(Array(piv)))
    Lv = [view(G, :, :, p) for p in 1:P]
    apply_S! = (dst, src) -> begin
        CUDA.CUBLAS.gemm_strided_batched!('N','N', ComplexF64(1), B3, src, ComplexF64(0), permbuf)
        _si_gather_rows!(dst, permbuf, perm)
        Dv = [view(dst, :, :, p) for p in 1:P]
        CUDA.CUBLAS.trsm_batched!('L','L','N','U', ComplexF64(1), Lv, Dv)
        CUDA.CUBLAS.trsm_batched!('L','U','N','N', ComplexF64(1), Lv, Dv)
        return dst
    end
    for _ in 1:Q
        apply_S!(Y, X)
        _si_cholqr!(Y, gram, rinv, tmp)
        X, Y = Y, X
    end
    apply_S!(Y, X)
    T = gram
    CUDA.CUBLAS.gemm_strided_batched!('C','N', ComplexF64(1), X, Y, ComplexF64(0), T)
    Th = Array(T)
    CUDA.unsafe_free!(G)
    Threads.@threads for p in 1:P
        μs = eigvals(@view Th[:, :, p])
        append!(cands[p], (σ + 1 / μ for μ in μs if abs(μ) > mu_tol))
    end
    return
end

# ── Batched contour-integral (Beyn) moment accumulation (see TJLF/src/tjlf_batched_si.jl) ─────
# Computes A_k = Σ_j wts[j,k]·(A - z_j B)⁻¹ B V over the whole batch with the same batched
# getrf/gather/trsm stack as the SI solver (no Xgeev, concurrency-safe). The pencil batch is
# sharded so the device footprint (A3+B3+G ~ 3n²·Ps plus the n×L blocks) fits in free memory;
# the per-pencil SVD/eig post-processing stays on the host in TJLF.contour_eigvals_batch.
function _contour_moments_shard!(mom, As, Bs, V, nodes, wts, rng)
    n, L = size(V)
    Ps = length(rng)
    K2 = size(wts, 2)
    A3h = Array{ComplexF64}(undef, n, n, Ps); B3h = similar(A3h)
    @inbounds for (q, p) in enumerate(rng)
        A3h[:, :, q] .= As[p]; B3h[:, :, q] .= Bs[p]
    end
    A3 = CUDA.CuArray(A3h); B3 = CUDA.CuArray(B3h)
    G  = similar(A3)
    V3 = CUDA.CuArray(repeat(reshape(V, n, L, 1), 1, 1, Ps))
    W3 = CUDA.zeros(ComplexF64, n, L, Ps)                       # B·V, node-independent
    Y  = CUDA.zeros(ComplexF64, n, L, Ps)
    Mk = [CUDA.zeros(ComplexF64, n, L, Ps) for _ in 1:K2]
    CUDA.CUBLAS.gemm_strided_batched!('N', 'N', ComplexF64(1), B3, V3, ComplexF64(0), W3)
    for (j, z) in enumerate(nodes)
        G .= A3 .- z .* B3
        piv, _, _ = CUDA.CUBLAS.getrf_strided_batched!(G, true)  # LU in place
        perm = CuArray(_si_ipiv_to_perm(Array(piv)))
        _si_gather_rows!(Y, W3, perm)
        Lv = [view(G, :, :, q) for q in 1:Ps]
        Dv = [view(Y, :, :, q) for q in 1:Ps]
        CUDA.CUBLAS.trsm_batched!('L', 'L', 'N', 'U', ComplexF64(1), Lv, Dv)
        CUDA.CUBLAS.trsm_batched!('L', 'U', 'N', 'N', ComplexF64(1), Lv, Dv)
        for k in 1:K2
            Mk[k] .+= wts[j, k] .* Y
        end
    end
    for k in 1:K2
        mom[:, :, rng, k] = Array(Mk[k])
    end
    foreach(CUDA.unsafe_free!, (A3, B3, G, V3, W3, Y))
    foreach(CUDA.unsafe_free!, Mk)
    return
end

function _contour_moments(As, Bs, V, nodes, wts)
    P = length(As)
    n, L = size(V)
    K2 = size(wts, 2)
    mom = Array{ComplexF64,4}(undef, n, L, P, K2)
    _with_device_slot() do
        bytes_per_p = (3 * n^2 + (3 + K2) * n * L) * sizeof(ComplexF64)
        Ps_max = clamp(floor(Int, 0.7 * CUDA.available_memory() / bytes_per_p), 1, P)
        for lo in 1:Ps_max:P
            _contour_moments_shard!(mom, As, Bs, V, nodes, wts, lo:min(lo + Ps_max - 1, P))
        end
    end
    return mom
end

# ── Coverage-certified adaptive SI: persistent batched handle (see tjlf_batched_si.jl) ────────
# One handle keeps the (A,B) stacks and work buffers on ONE pinned device across adaptive
# rounds; each solve(σs) is a batched shift-invert sweep with a PER-PENCIL shift vector
# (G = A - σ_p B broadcasts over the batch), returning per-pencil Ritz values AND their
# relative residuals in the ORIGINAL pencil, computed on device (only M×P norms cross PCIe).
function _csi_open(As, Bs, M::Int, Q::Int)
    P = length(As); n = size(As[1], 1)
    A3h = Array{ComplexF64}(undef, n, n, P); B3h = similar(A3h)
    @inbounds for p in 1:P
        A3h[:, :, p] .= As[p]; B3h[:, :, p] .= Bs[p]
    end
    devs = _devices()
    dev = devs[mod1(Threads.atomic_add!(_GPU_RR, 1) + 1, length(devs))]
    CUDA.device!(dev)
    id = CUDA.deviceid(dev)
    _warmup_device!(id)
    sem = _device_sem(id)

    A3 = CUDA.CuArray(A3h); B3 = CUDA.CuArray(B3h)
    bufs = (X = CUDA.zeros(ComplexF64, n, M, P), Y = CUDA.zeros(ComplexF64, n, M, P),
            permbuf = CUDA.zeros(ComplexF64, n, M, P), tmp = CUDA.zeros(ComplexF64, n, M, P),
            gram = CUDA.zeros(ComplexF64, M, M, P), rinv = CUDA.zeros(ComplexF64, M, M, P),
            U = CUDA.zeros(ComplexF64, M, M, P), Λ = CUDA.zeros(ComplexF64, M, M, P))
    rng = Random.MersenneTwister(0x51c0de)
    X0h = Matrix(qr(randn(rng, ComplexF64, n, M)).Q)
    X0rep = CUDA.CuArray(repeat(reshape(X0h, n, M, 1), 1, 1, P))

    solve = function (σs::Vector{ComplexF64})
        CUDA.device!(dev)
        Base.acquire(sem)
        try
            X, Y = bufs.X, bufs.Y
            copyto!(X, X0rep)
            σd = CUDA.CuArray(reshape(σs, 1, 1, P))
            G = A3 .- σd .* B3
            piv, _, _ = CUDA.CUBLAS.getrf_strided_batched!(G, true)
            perm = CuArray(_si_ipiv_to_perm(Array(piv)))
            Lv = [view(G, :, :, p) for p in 1:P]
            apply_S! = (dst, src) -> begin
                CUDA.CUBLAS.gemm_strided_batched!('N','N', ComplexF64(1), B3, src,
                                                  ComplexF64(0), bufs.permbuf)
                _si_gather_rows!(dst, bufs.permbuf, perm)
                Dv = [view(dst, :, :, p) for p in 1:P]
                CUDA.CUBLAS.trsm_batched!('L','L','N','U', ComplexF64(1), Lv, Dv)
                CUDA.CUBLAS.trsm_batched!('L','U','N','N', ComplexF64(1), Lv, Dv)
                return dst
            end
            for _ in 1:Q
                apply_S!(Y, X)
                _csi_cholqr!(Y, bufs.gram, bufs.rinv, bufs.tmp)
                X, Y = Y, X
            end
            apply_S!(Y, X)
            CUDA.CUBLAS.gemm_strided_batched!('C','N', ComplexF64(1), X, Y,
                                              ComplexF64(0), bufs.gram)
            Th = Array(bufs.gram)
            CUDA.unsafe_free!(G); CUDA.unsafe_free!(perm); CUDA.unsafe_free!(σd)

            # Host M×M eigs of T; degenerate lanes (NaN from a singular shift) are skipped.
            λh = fill(complex(NaN, NaN), M, P)
            Uh = zeros(ComplexF64, M, M, P); Λh = zeros(ComplexF64, M, M, P)
            Threads.@threads for p in 1:P
                E = try
                    eigen(@view Th[:, :, p])
                catch
                    nothing
                end
                E === nothing && continue
                Uh[:, :, p] .= E.vectors
                for i in 1:M
                    μ = E.values[i]
                    (isfinite(abs(μ)) && abs(μ) > 1e-14) || continue
                    λh[i, p] = σs[p] + 1 / μ
                    Λh[i, i, p] = λh[i, p]
                end
            end
            copyto!(bufs.U, Uh); copyto!(bufs.Λ, Λh)
            # Ritz residuals on device: x_i = X u_i, r_i = ‖Ax_i − λ_i Bx_i‖/(‖Ax_i‖+‖λ_i Bx_i‖).
            CUDA.CUBLAS.gemm_strided_batched!('N','N', ComplexF64(1), X, bufs.U,
                                              ComplexF64(0), Y)            # Y = XU
            CUDA.CUBLAS.gemm_strided_batched!('N','N', ComplexF64(1), A3, Y,
                                              ComplexF64(0), bufs.permbuf) # PA
            CUDA.CUBLAS.gemm_strided_batched!('N','N', ComplexF64(1), B3, Y,
                                              ComplexF64(0), bufs.tmp)     # PB
            CUDA.CUBLAS.gemm_strided_batched!('N','N', ComplexF64(1), bufs.tmp, bufs.Λ,
                                              ComplexF64(0), X)            # PBΛ
            na = sqrt.(Array(dropdims(sum(abs2, bufs.permbuf; dims = 1); dims = 1)))
            nb = sqrt.(Array(dropdims(sum(abs2, X; dims = 1); dims = 1)))
            bufs.permbuf .-= X
            nd = sqrt.(Array(dropdims(sum(abs2, bufs.permbuf; dims = 1); dims = 1)))
            resh = nd ./ max.(na .+ nb, eps())
            return λh, resh
        finally
            Base.release(sem)
        end
    end
    closed = Ref(false)
    close = () -> begin
        closed[] && return                    # idempotent: compaction closes, `finally` re-closes
        closed[] = true
        CUDA.device!(dev)
        foreach(CUDA.unsafe_free!, (A3, B3, X0rep))
        foreach(CUDA.unsafe_free!, values(bufs))
        return
    end
    return (solve = solve, close = close)
end

# Cholesky-QR with per-pencil failure isolation: a degenerate lane (NaN block from a singular
# shift) keeps its junk instead of crashing the whole batch; its residuals come out Inf.
function _csi_cholqr!(X, gram, rinv, tmp)
    _, mm, P = size(X)
    CUDA.CUBLAS.gemm_strided_batched!('C','N', ComplexF64(1), X, X, ComplexF64(0), gram)
    Gh = Array(gram); Rinv = similar(Gh)
    Threads.@threads for p in 1:P
        @views begin
            R = try
                cholesky(Hermitian((Gh[:, :, p] .+ Gh[:, :, p]') ./ 2)).U
            catch
                nothing
            end
            Rinv[:, :, p] .= R === nothing ? Matrix{ComplexF64}(I, mm, mm) : inv(R)
        end
    end
    copyto!(rinv, Rinv)
    CUDA.CUBLAS.gemm_strided_batched!('N','N', ComplexF64(1), X, rinv, ComplexF64(0), tmp)
    copyto!(X, tmp)
    return X
end

# Largest per-device batch for the certified-SI stack: A3+B3+G (3n²) + 4 n×M blocks + small M².
_csi_maxp(n::Int, M::Int) = begin
    mem = CUDA.totalmem(first(_devices()))
    max(1, floor(Int, 0.55 * mem / ((3 * n^2 + 4 * n * M + 4 * M^2) * sizeof(ComplexF64))))
end

function _si_batch(As, Bs, cfg)
    P = length(As); n = size(As[1], 1); M = cfg.M; Q = cfg.Q
    out = Vector{Vector{ComplexF64}}(undef, P)
    _with_device_slot() do
        # Memory-aware sharding over the pencil batch. Peak device footprint per pencil is
        # A3 + B3 + the per-shift LU buffer G (3 n²) plus the subspace/Gram buffers (X,Y,permbuf,tmp
        # = 4 nM; gram,rinv = 2 M²), all complex. A full ~256-pencil round at large n (UCP n=2400 →
        # ~0.28 GB/pencil) far exceeds a 40 GB A100, so cap the shard to ~60% of free memory and free
        # the device buffers between shards. At small n (DIII-D n≤1440) this is a single shard, so the
        # bit-for-bit behavior there is unchanged.
        bytes_per_p = (3 * n^2 + 4 * n * M + 2 * M^2) * sizeof(ComplexF64)
        Ps_max = clamp(floor(Int, 0.6 * CUDA.available_memory() / bytes_per_p), 1, P)
        X0h = Matrix(qr(randn(ComplexF64, n, M)).Q)
        for lo in 1:Ps_max:P
            rng = lo:min(lo + Ps_max - 1, P); Ps = length(rng)
            A3h = Array{ComplexF64}(undef, n, n, Ps); B3h = similar(A3h)
            @inbounds for (q, p) in enumerate(rng)
                A3h[:, :, q] .= As[p]; B3h[:, :, q] .= Bs[p]
            end
            A3 = CUDA.CuArray(A3h); B3 = CUDA.CuArray(B3h)
            bufs = (CUDA.zeros(ComplexF64, n, M, Ps), CUDA.zeros(ComplexF64, n, M, Ps),
                    CUDA.zeros(ComplexF64, n, M, Ps), CUDA.zeros(ComplexF64, n, M, Ps),
                    CUDA.zeros(ComplexF64, M, M, Ps), CUDA.zeros(ComplexF64, M, M, Ps))
            cands = [ComplexF64[] for _ in 1:Ps]
            for σ in cfg.shifts
                _si_shift!(cands, A3, B3, σ, X0h, bufs, M, Q, cfg.mu_tol)
            end
            for (q, p) in enumerate(rng)
                out[p] = TJLF._dedup_ritz(cands[q], cfg.dedup_tol)
            end
            CUDA.unsafe_free!(A3); CUDA.unsafe_free!(B3); foreach(CUDA.unsafe_free!, bufs)
        end
        return out
    end
end

function __init__()
    TJLF._CUDA_SI_BATCH[] = (As, Bs, cfg) -> _si_batch(As, Bs, cfg)
    TJLF._CUDA_CONTOUR_MOMENTS[] = (As, Bs, V, nodes, wts) -> _contour_moments(As, Bs, V, nodes, wts)
    TJLF._CUDA_CSI_OPEN[] = (As, Bs, M, Q) -> _csi_open(As, Bs, M, Q)
    TJLF._CUDA_CSI_MAXP[] = (n, M) -> _csi_maxp(n, M)
    TJLF._CUDA_FUNCTIONAL[]   = () -> CUDA.functional()
    TJLF._CUDA_DEVICE_COUNT[] = () -> length(CUDA.devices())
    # Returns Vector{ComplexF64} of eigenvalues of B⁻¹A, computed entirely on GPU
    # (CUSOLVER getrf/getrs for the solve, then Xgeev for the diagonalisation).
    # B is only needed for the LU factorisation; its overwritten form is unused.
    TJLF._CUDA_SOLVE[] = (A::Matrix{ComplexF64}, B::Matrix{ComplexF64}) -> begin
        _with_device_slot() do
            A_gpu = CUDA.CuArray(A)
            B_gpu = CUDA.CuArray(B)
            ipiv  = CUDA.CuArray{Int32}(undef, size(B_gpu, 1))
            B_factored, ipiv, _ = CUDA.CUSOLVER.getrf!(B_gpu, ipiv)
            A_gpu = CUDA.CUSOLVER.getrs!('N', B_factored, ipiv, A_gpu)
            W, _, _ = CUDA.CUSOLVER.Xgeev!('N', 'N', A_gpu)
            return Array(W)
        end
    end

    # Eigenvector inverse-iteration step on GPU: solve Z x = b (LU via getrf/getrs).
    # Mirrors _CUDA_SOLVE's device selection + per-device locking so it composes with the
    # eigenvalue solve (same device, serialized intra-process, overlapped across MPS clients).
    TJLF._CUDA_LU_SOLVE[] = (Z::Matrix{ComplexF64}, b::Vector{ComplexF64}) -> begin
        _with_device_slot() do
            Z_gpu = CUDA.CuArray(Z)
            b_gpu = CUDA.CuArray(reshape(b, :, 1))
            ipiv  = CUDA.CuArray{Int32}(undef, size(Z_gpu, 1))
            Z_factored, ipiv, _ = CUDA.CUSOLVER.getrf!(Z_gpu, ipiv)
            x_gpu = CUDA.CUSOLVER.getrs!('N', Z_factored, ipiv, b_gpu)
            return vec(Array(x_gpu))
        end
    end

    # Complex{Dual} eigensolve: eigenvalues of M=B⁻¹A plus IFT derivatives ∂λᵢ for each
    # partial. Mirrors the CPU Dual kernel (TJLFForwardDiffExt) but keeps every O(n³) step
    # on the GPU: getrf/getrs for M and ∂M (one LU of B reused), Xgeev('N','V') for the
    # eigenpairs, then ∂λᵢ=(R⁻¹∂M R)[i,i]. Round-robin device + per-device lock match
    # _CUDA_SOLVE so team workers overlap across GPUs via MPS while staying serialized
    # within a CUDA context (Xgeev is not concurrency-safe intra-context).
    TJLF._CUDA_SOLVE_EIG_GRAD[] = (A::Matrix{ComplexF64}, B::Matrix{ComplexF64},
                                   dA::Vector{Matrix{ComplexF64}}, dB::Vector{Matrix{ComplexF64}}) -> begin
        _with_device_slot() do
            np = length(dA)
            m  = size(A, 1)
            A_gpu = CUDA.CuArray(A)
            B_gpu = CUDA.CuArray(B)
            ipivB = CUDA.CuArray{Int32}(undef, m)
            B_fact, ipivB, _ = CUDA.CUSOLVER.getrf!(B_gpu, ipivB)
            # Mf = B⁻¹A (getrs! overwrites A_gpu and returns it).
            Mf = CUDA.CUSOLVER.getrs!('N', B_fact, ipivB, A_gpu)
            # ∂M[k] = B⁻¹(∂A[k] − ∂B[k]·Mf), reusing the LU of B. Done BEFORE Xgeev overwrites Mf.
            dMs = Vector{CUDA.CuMatrix{ComplexF64}}(undef, np)
            for k in 1:np
                dAk = CUDA.CuArray(dA[k])
                dBk = CUDA.CuArray(dB[k])
                rhs = dAk - dBk * Mf
                dMs[k] = CUDA.CUSOLVER.getrs!('N', B_fact, ipivB, rhs)
            end
            # Eigenpairs of Mf: eigenvalues W + right vectors VR (Xgeev! overwrites Mf).
            W, _, VR = CUDA.CUSOLVER.Xgeev!('N', 'V', Mf)
            # Factor a COPY of VR (getrf! is in-place); keep VR itself for the ∂M·R product.
            VRf = copy(VR)
            ipivR = CUDA.CuArray{Int32}(undef, m)
            VRf, ipivR, _ = CUDA.CUSOLVER.getrf!(VRf, ipivR)
            dλ_re = zeros(Float64, m, np)
            dλ_im = zeros(Float64, m, np)
            for k in 1:np
                C  = CUDA.CUSOLVER.getrs!('N', VRf, ipivR, dMs[k] * VR)  # = R⁻¹ ∂M R
                Ch = Array(C)
                @inbounds for i in 1:m
                    dλ_re[i, k] = real(Ch[i, i])
                    dλ_im[i, k] = imag(Ch[i, i])
                end
            end
            return (Array(W), dλ_re, dλ_im)
        end
    end
end

end

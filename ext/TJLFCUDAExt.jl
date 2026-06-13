module TJLFCUDAExt

using CUDA
using TJLF

# Per-device dispatch locks. cusolverDnCreate races when many threads call it
# simultaneously on first use, so Julia-level dispatch is serialized — but only
# *per device*. Using one lock per GPU (rather than a single global lock) lets
# solves on different GPUs proceed concurrently, which is what enables
# multi-GPU-per-radius scaling: e.g. with --gpus-per-task=2 the round-robin below
# splits a radius's ~1024 eigensolves across both A100s, ~halving per-radius time.
# Within a single device, dispatch stays serialized (the GPU is internally
# parallel, so this has negligible cost). With one visible GPU this reduces to a
# single lock — identical to the previous behavior.
const _LOCKS_GUARD = ReentrantLock()
const _GPU_DEVICE_LOCKS = Dict{Int,ReentrantLock}()
function _device_lock(id::Integer)
    lock(_LOCKS_GUARD) do
        get!(() -> ReentrantLock(), _GPU_DEVICE_LOCKS, Int(id))
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
    lock(_LOCKS_GUARD) do
        devs = _DEVICES[]
        isempty(devs) || return devs
        devs = collect(CUDA.devices())
        _DEVICES[] = devs
        return devs
    end
end

function __init__()
    TJLF._CUDA_FUNCTIONAL[]   = () -> CUDA.functional()
    TJLF._CUDA_DEVICE_COUNT[] = () -> length(CUDA.devices())
    # Returns Vector{ComplexF64} of eigenvalues of B⁻¹A, computed entirely on GPU
    # (CUSOLVER getrf/getrs for the solve, then Xgeev for the diagonalisation).
    # B is only needed for the LU factorisation; its overwritten form is unused.
    TJLF._CUDA_SOLVE[] = (A::Matrix{ComplexF64}, B::Matrix{ComplexF64}) -> begin
        devs = _devices()
        n = length(devs)
        # round-robin device selection across all visible GPUs
        dev = devs[mod1(Threads.atomic_add!(_GPU_RR, 1) + 1, n)]
        CUDA.device!(dev)
        lock(_device_lock(CUDA.deviceid(dev))) do
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
        devs = _devices()
        n = length(devs)
        dev = devs[mod1(Threads.atomic_add!(_GPU_RR, 1) + 1, n)]
        CUDA.device!(dev)
        lock(_device_lock(CUDA.deviceid(dev))) do
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
        devs = _devices()
        ndev = length(devs)
        dev = devs[mod1(Threads.atomic_add!(_GPU_RR, 1) + 1, ndev)]
        CUDA.device!(dev)
        lock(_device_lock(CUDA.deviceid(dev))) do
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

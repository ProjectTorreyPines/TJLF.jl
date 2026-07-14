# Batched shift-invert dominant-mode eigensolver.
#
# Motivation (measured on real TJLF pencils, see TJLFEP/build/ad/benchmark_*_eigs*.jl):
#   * cuSOLVER Xgeev is ~98% of every GPU eigensolve, cannot be batched, and is NOT
#     concurrency-safe within a CUDA context — so the dense GPU path is forced serial.
#   * TJLF only consumes the few MOST-UNSTABLE modes (tjlf_LINEAR_SOLUTION sorts by growth
#     rate and keeps NMODES / the per-branch leaders). Those live in a small cluster near the
#     origin (gamma in (0, ~0.2], |freq| <~ 1.5) while the spectrum bulk is strongly damped.
#   * A shift-invert subspace iteration targeting that window is fully expressible as BATCHED
#     cuBLAS primitives (batched LU + triangular solves + GEMM + tiny host M×M eig), which IS
#     concurrency-safe and fills the GPU: measured ~6x faster per pencil than serial Xgeev at
#     n=720, with the physical branch leaders recovered to ~1e-6.
#
# This is an APPROXIMATE eigensolver (dominant modes only). Use it for eigenvalue-only stages
# (IFLUX=false onset scans); confirm the winning configuration with the dense solver.

"""
    SIConfig(; shifts, M=16, Q=12, dedup_tol=1e-8, mu_tol=1e-10)

Configuration for the batched shift-invert solver. `shifts` are complex targets σ placed to
cover the physically relevant unstable window; modes near each σ become dominant under
S = (A - σB)⁻¹B and are found by `Q` steps of `M`-dimensional subspace iteration, then
back-transformed λ = σ + 1/μ and unioned/deduped across shifts.
"""
Base.@kwdef struct SIConfig
    # Shifts cover the two physical branches (ion: imag(lambda) up to ~1.5; electron: a dense
    # cluster hugging the real axis, |freq| <~ 0.12). The near-axis band (+/-0.05, +/-0.12) is
    # sampled densely because both weak ion modes and the entire electron branch live there — a
    # sparse jump straight to +/-0.25 was missing ~5/37 electron leaders at n=1440.
    shifts::Vector{ComplexF64} = ComplexF64[0.02,
                                            0.02+0.05im, 0.02-0.05im,
                                            0.02+0.12im, 0.02-0.12im,
                                            0.02+0.25im, 0.02-0.25im,
                                            0.05+0.6im,  0.05-0.6im,
                                            0.05+1.1im,  0.05-1.1im,
                                            0.05+1.5im,  0.05-1.5im]
    M::Int = 16
    Q::Int = 12
    dedup_tol::Float64 = 1e-8
    mu_tol::Float64 = 1e-10
end

# Union/dedup Ritz values (sorted by (re, im)), collapsing near-duplicates found by multiple shifts.
function _dedup_ritz(lams::AbstractVector{ComplexF64}, tol::Float64)
    isempty(lams) && return ComplexF64[]
    s = sort(lams; by = x -> (real(x), imag(x)))
    out = ComplexF64[s[1]]
    @inbounds for l in @view s[2:end]
        (abs(l - out[end]) > tol * max(1.0, abs(l))) && push!(out, l)
    end
    return out
end

# CPU reference: per-pencil shift-invert Arnoldi (KrylovKit) unioned over shifts. Correctness/
# fallback path — NOT fast (sequential dense LU per shift), used for testing and when no GPU.
function _si_eigvals_batch_cpu(As::AbstractVector{<:AbstractMatrix{ComplexF64}},
                               Bs::AbstractVector{<:AbstractMatrix{ComplexF64}}, cfg::SIConfig)
    P = length(As)
    out = Vector{Vector{ComplexF64}}(undef, P)
    for p in 1:P
        A = As[p]; B = Bs[p]; n = size(A, 1)
        lams = ComplexF64[]
        for σ in cfg.shifts
            F = lu(A - σ * B)
            op = x -> F \ (B * x)
            x0 = ones(ComplexF64, n) ./ sqrt(n)
            μs, _, _ = KrylovKit.eigsolve(op, x0, cfg.M, :LM;
                                          krylovdim = max(2cfg.M, 30), tol = 1e-9,
                                          maxiter = 30, ishermitian = false, verbosity = 0)
            append!(lams, (σ + 1 / μ for μ in μs if abs(μ) > cfg.mu_tol))
        end
        out[p] = _dedup_ritz(lams, cfg.dedup_tol)
    end
    return out
end

"""
    si_eigvals_batch(As, Bs; cfg=SIConfig(), use_gpu=false) -> Vector{Vector{ComplexF64}}

Approximate dominant eigenvalues of each generalized pencil `(As[p], Bs[p])` (standard problem
M = B⁻¹A) near the configured shifts. Returns, per pencil, the deduped back-transformed Ritz
values (a superset of the unstable modes). On `use_gpu=true` with the CUDA extension loaded this
runs as one batched sweep on the GPU; otherwise it falls back to the CPU reference.
"""
function si_eigvals_batch(As::AbstractVector{<:AbstractMatrix{ComplexF64}},
                          Bs::AbstractVector{<:AbstractMatrix{ComplexF64}};
                          cfg::SIConfig = SIConfig(), use_gpu::Bool = false)
    length(As) == length(Bs) || throw(DimensionMismatch("As and Bs must have equal length"))
    isempty(As) && return Vector{Vector{ComplexF64}}()
    if use_gpu && _CUDA_SI_BATCH[] !== nothing
        return _CUDA_SI_BATCH[](As, Bs, cfg)
    end
    return _si_eigvals_batch_cpu(As, Bs, cfg)
end

# ══ Contour-integral (Beyn) batched eigensolver ═══════════════════════════════════════════════
#
# Second-chance replacement for the point-shift SI above. The fixed-shift scheme misses leaders
# whose eigenvalues sit where no shift was placed (IR101: 333/503 ion leaders), and the
# geev-calibrated adaptive scheme still misses modes absent from the calibration sample — both
# are COVERAGE-BY-SAMPLING failures. A contour integral around the physically relevant window
# finds ALL eigenvalues inside the region: coverage is geometric, not sampled.
#
# Method (Sakurai–Sugiura Rayleigh–Ritz, SS-RR): for the linear pencil T(z) = A - zB, the
# filtered moment blocks
#     A_k = (1/2πi) ∮_Γ φ(z)^k T(z)⁻¹ B V dz ,   φ(z) = (z - c)/ρ ,  V ∈ C^{n×L} random
# span (up to quadrature leakage) the invariant subspace of the eigenvalues INSIDE Γ. We
# orthogonalize the stacked moments [A_0 … A_{2K-1}] with a rank-revealing SVD and Rayleigh-Ritz
# the ORIGINAL pencil (UᴴAU, UᴴBU) — unlike the Hankel variant (Beyn Alg. 2) this stays accurate
# when weakly damped modes crowd the contour (real TJLF spectra have ~100 modes hugging γ≲0),
# because quadrature leakage only wastes subspace capacity instead of corrupting the small
# eigenproblem. Ritz pairs are then residual-CERTIFIED against the original (A,B), and
# near-boundary candidates are polished by Rayleigh-quotient iteration.
#
# Every expensive step is a batched cuBLAS primitive over all P pencils (same getrf/trsm/gemm
# stack as the SI solver — no Xgeev, concurrency-safe); the SVD + projected eig are per-pencil
# host work. Failure modes are DETECTABLE, not silent:
#   * rank saturation (m ≥ sat_frac·2K·L): subspace capacity may be exhausted → flag pencil;
#   * an un-certifiable candidate in the physically relevant half-window → flag pencil;
#   * a refined candidate that moved suspiciously far (identity uncertain) → flag pencil.
# Flagged pencils are routed to the exact dense eigensolve (per pencil, not per radius).

"""
    ContourConfig(; re_lo=-0.02, re_hi=0.8, im_max=2.6, n_long=18, n_short=6, L=64, K=3, ...)

Configuration for [`contour_eigvals_batch`](@ref). The contour is the RECTANGLE
`[re_lo, re_hi] × [-im_max, im_max]` in the complex λ plane — a rectangle (not an ellipse)
because the physically consumed AE/ion modes live in the corners (small growth rate, large
|frequency|), which an ellipse would cut off. `re_lo` must sit below 0 so marginal (γ≈0)
modes are interior; `re_hi`/`im_max` bound the strongest expected drive and mode frequency
(audit against real spectra — a consumed mode outside the window is invisible BY DESIGN).

- `n_long`/`n_short`: Gauss-Legendre nodes on each long (vertical, length `2·im_max`) / short
  (horizontal) side; total factorizations per batch = `2(n_long + n_short)`.
- `L`: probe block columns; `K`: block moments. Subspace capacity = `2K·L` must exceed the
  in-window eigenvalue count INCLUDING the weakly damped modes crowding just inside `re_lo`
  (measured real TJLF pencils at n=720: up to ~110 for `re_lo=-0.02`) AND the quadrature
  leakage of near-boundary exterior modes.
- `rank_tol`: relative SVD truncation of the block-Hankel moment matrix.
- `resid_tol`: relative residual `‖Ax-λBx‖/(‖Ax‖+|λ|‖Bx‖)` below which a candidate is certified.
- `refine`: max Rayleigh-quotient iterations to polish uncertified candidates with
  `real(λ) > cert_re` (quadrature only ENUMERATES near-boundary modes coarsely; refinement
  makes them precise, certification makes them trusted).
- `cert_re`: only candidates with `real(λ) > cert_re` must certify — the deep-stable crowd is
  physically inert (never a branch leader, never kept) and is returned unrefined.
- `max_move`: a refined candidate moving farther than this may have jumped to a DIFFERENT
  eigenvalue (the original one would then be silently lost) → flag the pencil.
- `sat_frac`: rank ≥ `sat_frac·2K·L` flags the pencil (subspace capacity possibly saturated).
"""
Base.@kwdef struct ContourConfig
    re_lo::Float64 = -0.02
    re_hi::Float64 = 0.8
    im_max::Float64 = 2.6
    n_long::Int = 18
    n_short::Int = 6
    L::Int = 64
    K::Int = 3
    rank_tol::Float64 = 1e-8
    resid_tol::Float64 = 1e-7
    refine::Int = 2
    cert_re::Float64 = -0.005
    max_move::Float64 = 0.05
    sat_frac::Float64 = 0.9
end

_contour_inside(λ::ComplexF64, cfg::ContourConfig) =
    cfg.re_lo ≤ real(λ) ≤ cfg.re_hi && abs(imag(λ)) ≤ cfg.im_max

# Reference point / scale of the window, used for the moment monomials φ(z) = (z-c)/ρ.
_contour_scale(cfg::ContourConfig) =
    (complex((cfg.re_lo + cfg.re_hi) / 2, 0.0), max((cfg.re_hi - cfg.re_lo) / 2, cfg.im_max))

# Quadrature nodes z_j on the rectangle (counterclockwise, composite Gauss-Legendre per side,
# panels of ≤ 8 nodes) and weights w[j,k] such that A_{k-1} ≈ Σ_j w[j,k] · T(z_j)⁻¹ B V,
# approximating A_{k-1} = (1/2πi) ∮ φ(z)^{k-1} T(z)⁻¹ B V dz.
function _contour_nodes(cfg::ContourConfig)
    c, ρ = _contour_scale(cfg)
    corners = ComplexF64[complex(cfg.re_hi, -cfg.im_max), complex(cfg.re_hi, cfg.im_max),
                         complex(cfg.re_lo, cfg.im_max), complex(cfg.re_lo, -cfg.im_max)]
    z = ComplexF64[]
    dzw = ComplexF64[]                                # GL weight × dz/dt (t ∈ [-1,1] per panel)
    for s in 1:4
        za, zb = corners[s], corners[mod1(s + 1, 4)]
        npts = iseven(s) ? cfg.n_short : cfg.n_long   # sides 1,3 vertical (long); 2,4 horizontal
        npanel = cld(npts, 8)
        g = cld(npts, npanel)
        x, w = FastGaussQuadrature.gausslegendre(g)
        for q in 1:npanel
            pa = za + (q - 1) / npanel * (zb - za)
            pb = za + q / npanel * (zb - za)
            half = (pb - pa) / 2
            for i in 1:g
                push!(z, (pa + pb) / 2 + half * x[i])
                push!(dzw, w[i] * half)
            end
        end
    end
    N = length(z)
    wts = Matrix{ComplexF64}(undef, N, 2 * cfg.K)
    for j in 1:N
        φ = (z[j] - c) / ρ
        for k in 1:2*cfg.K
            wts[j, k] = φ^(k - 1) * dzw[j] / (2π * im)
        end
    end
    return z, wts
end

# CPU reference for the moment accumulation (test/fallback path; the GPU hook computes the same
# quantity with batched getrf/trsm/gemm). Returns host (n, L, P, 2K).
function _contour_moments_cpu(As, Bs, V::Matrix{ComplexF64},
                              nodes::Vector{ComplexF64}, wts::Matrix{ComplexF64})
    P = length(As); n, L = size(V); K2 = size(wts, 2)
    mom = zeros(ComplexF64, n, L, P, K2)
    Threads.@threads for p in 1:P
        W = Bs[p] * V
        for (j, z) in enumerate(nodes)
            Y = lu!(As[p] - z * Bs[p]) \ W
            for k in 1:K2
                @views mom[:, :, p, k] .+= wts[j, k] .* Y
            end
        end
    end
    return mom
end

# Per-pencil SS-RR extraction: orthonormalize the stacked moments [A_0 … A_{2K-1}] (n × 2KL)
# with a rank-revealing SVD, then Rayleigh-Ritz the ORIGINAL pencil on that subspace
# (Ã = UᴴAU, B̃ = UᴴBU, generalized eig). Returns (Ritz values, Ritz vectors, rank, saturated).
function _contour_extract(A::AbstractMatrix{ComplexF64}, B::AbstractMatrix{ComplexF64},
                          mom_p::AbstractArray{ComplexF64,3}, cfg::ContourConfig)
    n, L, K2 = size(mom_p)
    S = svd!(reshape(Array(mom_p), n, L * K2))
    m = count(>(cfg.rank_tol * max(S.S[1], eps())), S.S)
    m == 0 && return ComplexF64[], Matrix{ComplexF64}(undef, n, 0), 0, false
    sat = m ≥ cfg.sat_frac * L * K2
    U = S.U[:, 1:m]
    F = eigen!(U' * (A * U), U' * (B * U))
    ok = findall(isfinite, F.values)          # B̃ near-singular directions yield Inf — drop
    return F.values[ok], U * F.vectors[:, ok], m, sat
end

_contour_resid(A, B, λ::ComplexF64, x::AbstractVector{ComplexF64}) = begin
    Ax = A * x; Bx = B * x
    norm(Ax .- λ .* Bx) / max(norm(Ax) + abs(λ) * norm(Bx), eps())
end

# Rayleigh-quotient / inverse iteration polish of one candidate (host, one n×n LU per step).
function _contour_refine(A, B, λ::ComplexF64, x::Vector{ComplexF64}, cfg::ContourConfig)
    for _ in 1:cfg.refine
        w = try
            lu!(A - λ * B) \ (B * x)
        catch                                  # exactly singular shift: λ is (numerically) exact
            break
        end
        x = w ./ norm(w)
        λ = (x' * (A * x)) / (x' * (B * x))
    end
    return λ, x
end

"""
    contour_eigvals_batch(As, Bs; cfg=ContourConfig(), use_gpu=false, dense_fallback=:cpu)
        -> (vals::Vector{Vector{ComplexF64}}, flagged::BitVector, ranks::Vector{Int})

Eigenvalues of each generalized pencil `(As[p], Bs[p])` inside the configured rectangular
window, via a Sakurai–Sugiura contour-integral subspace + Rayleigh-Ritz on the original pencil.
All pencils must share one size `n`. Every returned eigenvalue is residual-certified
(`‖Ax-λBx‖/(‖Ax‖+|λ|‖Bx‖) ≤ resid_tol`; uncertified candidates in the physically relevant
half-window `real(λ) > cert_re` are Rayleigh-quotient refined first). Pencils whose failure
indicators fire (rank saturation, an uncertifiable relevant candidate, a refined candidate of
uncertain identity) are `flagged` and — unless `dense_fallback === :none` — solved exactly with
the dense eigensolver (`:cpu` threaded LAPACK, or `:gpu` serial cuSOLVER Xgeev), in which case
`vals[p]` is the full dense spectrum for those pencils.
"""
function contour_eigvals_batch(As::AbstractVector{<:AbstractMatrix{ComplexF64}},
                               Bs::AbstractVector{<:AbstractMatrix{ComplexF64}};
                               cfg::ContourConfig = ContourConfig(), use_gpu::Bool = false,
                               dense_fallback::Symbol = :cpu)
    length(As) == length(Bs) || throw(DimensionMismatch("As and Bs must have equal length"))
    P = length(As)
    P == 0 && return Vector{Vector{ComplexF64}}(), falses(0), Int[]
    n = size(As[1], 1)
    all(A -> size(A, 1) == n, As) || throw(DimensionMismatch("all pencils must share one size n"))

    # Deterministic orthonormal probe block: same V on every call for bit-reproducible scans.
    rng = Random.MersenneTwister(0x7ea1c0de)
    V = Matrix(qr(randn(rng, ComplexF64, n, cfg.L)).Q)
    nodes, wts = _contour_nodes(cfg)

    mom = (use_gpu && _CUDA_CONTOUR_MOMENTS[] !== nothing) ?
          _CUDA_CONTOUR_MOMENTS[](As, Bs, V, nodes, wts) :
          _contour_moments_cpu(As, Bs, V, nodes, wts)

    vals = Vector{Vector{ComplexF64}}(undef, P)
    flagged = falses(P)
    ranks = Vector{Int}(undef, P)
    Threads.@threads for p in 1:P
        A = As[p]; B = Bs[p]
        λs, X, m, sat = _contour_extract(A, B, view(mom, :, :, p, :), cfg)
        ranks[p] = m
        keep = ComplexF64[]
        bad = false
        for (i, λ0) in enumerate(λs)
            _contour_inside(λ0, cfg) || continue
            x = X[:, i]
            r = _contour_resid(A, B, λ0, x)
            λ = λ0
            if r > cfg.resid_tol
                # Deep-stable crowd (never consumed downstream): drop uncertified silently.
                real(λ0) > cfg.cert_re || continue
                λ, x = _contour_refine(A, B, λ0, x, cfg)
                r = _contour_resid(A, B, λ, x)
                if r > cfg.resid_tol || abs(λ - λ0) > cfg.max_move
                    bad = true                 # relevant candidate we could not pin down
                    continue
                end
            end
            _contour_inside(λ, cfg) && push!(keep, λ)
        end
        flagged[p] = sat || bad
        vals[p] = _dedup_ritz(keep, 1e-10)
    end

    _dense_fallback!(vals, flagged, As, Bs, dense_fallback)
    return vals, flagged, ranks
end

# Replace vals[p] with the exact dense spectrum for every flagged pencil.
function _dense_fallback!(vals, flagged, As, Bs, dense_fallback::Symbol)
    dense_fallback === :none && return vals
    idx = findall(flagged)
    isempty(idx) && return vals
    if dense_fallback === :gpu && _CUDA_SOLVE[] !== nothing
        for p in idx      # serial: Xgeev is not concurrency-safe intra-context
            vals[p] = _gpu_solve!(copy(As[p]), copy(Bs[p]))
        end
    else
        Threads.@threads for p in idx
            A2, B2 = copy(As[p]), copy(Bs[p])
            (A3, _, _) = gesv!(B2, A2)
            vals[p] = geev!('N', 'N', A3)[1]
        end
    end
    return vals
end

# ══ Coverage-certified adaptive multi-shift SI ════════════════════════════════════════════════
#
# Measured on real TJLF pencils, the contour (SS-RR) solver above is SAFE but ineffective: the
# spectrum has a quasi-continuum of ~100+ weakly damped modes hugging γ≲0, so any contour
# boundary near the axis runs through a mode crowd — quadrature leakage saturates the subspace
# (ranks ~350/384 at 48 nodes) and almost nothing certifies. There is no spectral gap to place a
# contour in. Shift-invert does not need a gap: it resolves the modes NEAREST each shift
# regardless of crowding. Its historical failure was coverage-by-sampling (fixed shifts missed
# 333/503 IR101 ion leaders; geev-calibrated adaptive shifts still missed ~8/1024).
#
# This solver keeps the fast batched SI kernel but makes coverage SELF-CERTIFYING, per pencil:
#   1. residual-certify every Ritz pair against the original pencil (batched on device):
#      only λ with ‖Ax-λBx‖/(‖Ax‖+|λ|‖Bx‖) ≤ resid_tol are returned or trusted;
#   2. each shift σ certifiably "owns" the disk D(σ, r), r = distance to its trust-th certified
#      Ritz value (subspace iteration converges the dominant = nearest modes first, so the
#      nearest `trust ≤ M/2` certified values are, with overwhelming probability, THE nearest
#      true modes — capped additionally at the nearest uncertified Ritz value);
#   3. the unstable window [0,re_hi]×[±im_max] (all downstream-consumed modes have γ>0) must be
#      covered by the union of disks; uncovered pencils get a NEW shift at their most uncovered
#      point next round — per-pencil shifts, still ONE batched sweep for the whole batch
#      (G = A - σ_p B broadcasts a shift vector);
#   4. pencils still uncovered after max_rounds are flagged → exact dense fallback (per pencil).
# No calibration geev is needed; coverage adapts to each pencil's own spectrum.

"""
    CertifiedSIConfig(; re_hi=0.8, im_max=2.6, row_re=0.02, row_dy=0.5, M=24, Q=12, ...)

Configuration for [`certified_si_eigvals_batch`](@ref).

- `re_hi`, `im_max`: the unstable window `[0,re_hi] × [-im_max,im_max]` that must be certifiably
  covered (audit against real spectra: consumed modes outside it are invisible BY DESIGN).
- `row_re`, `row_dy`, `row_dense_w`: initial shift row at `Re = row_re` hugging the axis (the
  electron branch and weak ion modes live there), dense (`row_dy`) only inside the crowd band
  `|ω| ≤ row_dense_w`; wider windows get `row2_dy`-spaced wing shifts, with adaptive rounds
  (batch-compacted, so cheap) mopping up per-pencil leftovers.
- `row2_re`, `row2_dy`: sparse second base row through the mid-window (usually mode-free in
  real spectra: its emptiness disks bulk-cover the right half instead of spending one adaptive
  round per hole).
- `M`, `Q`: subspace block size / iterations per shift (per pencil per shift, the `M` nearest
  modes are converged; only ~`trust` of them are relied on for coverage).
- `resid_tol`: relative residual below which a Ritz VALUE is certified (trusted downstream).
- `enum_tol`: looser residual below which a Ritz pair ENUMERATES a real nearby mode. Coverage
  radii are built from enumerated modes (a Ritz pair at residual 1e-3 locates a mode to
  ~κ·1e-3 even though its value is not yet certified); with only the certified set the disks
  collapse in the axis-hugging crowd, where full 1e-8 convergence of many modes is slow.
- `trust`: coverage radius = distance from σ to its `trust`-th enumerated Ritz value (also
  capped at the nearest NON-enumerated Ritz value); `trust ≤ M/2` leaves convergence margin.
- `min_cert`: fewer enumerated values than this ⇒ the shift claims no enumeration disk.
- `empty_frac`: emptiness disk. Block power iteration converges the DOMINANT eigenvalue of the
  shift-invert operator quickly, so the nearest Ritz distance `d` estimates the distance to the
  nearest true mode even when nothing certifies (empty window regions with no |μ| separation);
  `D(σ, empty_frac·d)` is claimed mode-free and therefore covered. Non-normality errs safe
  here: Ritz values live in the field of values, which OVERestimates the spectral radius, so
  the disk shrinks, never grows.
- `refine`, `cert_lo`, `max_move`: enumerated-but-uncertified candidates with
  `real(λ) > cert_lo` (the downstream-relevant ones) are polished by `refine` inverse/Rayleigh
  iterations; if they still fail `resid_tol`, or move by more than `max_move` (identity
  uncertain), the pencil is flagged. Uncertified crowd modes (`real(λ) ≤ cert_lo`) are dropped
  silently — they are never consumed downstream.
- `max_refine`: if a pencil accumulates more than this many uncertified relevant candidates,
  flag it instead of refining — one dense `geev` is cheaper and exact than dozens of n³ LU
  refinements (bounds host cost on wide windows / crowded spectra).
- `max_rounds`: adaptive rounds after the initial row before uncovered pencils are flagged.
- `x_fine`, `dx_fine`, `dy_fine`, `dx_coarse`, `dy_coarse`: coverage test CELLS over the
  window. A disk only covers a cell if it covers the ENTIRE cell (center distance + half cell
  diagonal ≤ radius), so the certificate has no inter-point gaps; the γ ∈ [0, x_fine] strip is
  resolved at `dx_fine × dy_fine` because the axis-hugging mode crowd forces small disks there.
"""
Base.@kwdef struct CertifiedSIConfig
    re_hi::Float64 = 0.8
    im_max::Float64 = 2.6
    row_re::Float64 = 0.02
    row_dy::Float64 = 0.25
    row_dense_w::Float64 = 2.7
    row2_re::Float64 = 0.45
    row2_dy::Float64 = 1.0
    M::Int = 32
    Q::Int = 20
    resid_tol::Float64 = 1e-8
    enum_tol::Float64 = 1e-2
    trust::Int = 16
    min_cert::Int = 2
    empty_frac::Float64 = 0.5
    refine::Int = 5
    cert_lo::Float64 = -0.005
    max_move::Float64 = 0.05
    max_refine::Int = 12
    max_rounds::Int = 60
    x_fine::Float64 = 0.06
    dx_fine::Float64 = 0.006
    dy_fine::Float64 = 0.012
    dx_coarse::Float64 = 0.025
    dy_coarse::Float64 = 0.05
    mu_tol::Float64 = 1e-12
    dedup_tol::Float64 = 1e-8
end

# CPU reference for one batched shift sweep: for each pencil p, shift-invert subspace iteration
# at σs[p], returning Ritz values λ[i,p] and relative residuals res[i,p] (i = 1..M).
# Same math as the GPU handle in TJLFCUDAExt; used for tests and when no GPU is available.
function _csi_open_cpu(As, Bs, M::Int, Q::Int)
    P = length(As); n = size(As[1], 1)
    rng = Random.MersenneTwister(0x51c0de)
    X0 = Matrix(qr(randn(rng, ComplexF64, n, M)).Q)
    solve = function (σs::Vector{ComplexF64})
        λ = fill(complex(NaN, NaN), M, P)
        res = fill(Inf, M, P)
        Threads.@threads for p in 1:P
            A = As[p]; B = Bs[p]
            F = try
                lu(A - σs[p] * B)
            catch
                continue                       # singular shift: no results from this pencil
            end
            X = copy(X0)
            for _ in 1:Q
                X = Matrix(qr!(F \ (B * X)).Q)
            end
            Y = F \ (B * X)
            T = X' * Y
            E = eigen!(T)
            XU = X * E.vectors
            AXU = A * XU; BXU = B * XU
            for i in 1:M
                μ = E.values[i]
                abs(μ) > 1e-14 || continue
                λi = σs[p] + 1 / μ
                a = @view AXU[:, i]; b = @view BXU[:, i]
                r = norm(a .- λi .* b) / max(norm(a) + abs(λi) * norm(b), eps())
                λ[i, p] = λi; res[i, p] = r
            end
        end
        return λ, res
    end
    return (solve = solve, close = () -> nothing)
end

# Per-pencil coverage state: a paving of the window by CELLS (center + half-diagonal). A disk
# covers a cell only when it contains the whole cell, so certified coverage has no gaps between
# test points (a real IR pencil lost a crowd mode exactly that way with a point grid).
# Only UNCOVERED cells are stored (vectors compact as disks land), and each keeps its running
# min-distance to the shifts placed so far — the wide IM_MAX=7 window has ~20k cells/pencil, and
# a full cells×shifts rescan per adaptive round starved the GPUs from the host side.
mutable struct _CSICover
    centers::Vector{ComplexF64}   # uncovered cell centers
    hdiag::Vector{Float64}        # matching half cell diagonals
    dmin::Vector{Float64}         # matching min distance to all shifts seen so far
end

function _CSICover(cfg::CertifiedSIConfig)
    centers = ComplexF64[]
    hdiag = Float64[]
    function pave!(x0, x1, dx, dy)
        nx = max(1, ceil(Int, (x1 - x0) / dx)); ny = max(1, ceil(Int, 2 * cfg.im_max / dy))
        ddx = (x1 - x0) / nx; ddy = 2 * cfg.im_max / ny
        h = hypot(ddx, ddy) / 2
        for jx in 1:nx, jy in 1:ny
            push!(centers, complex(x0 + (jx - 0.5) * ddx, -cfg.im_max + (jy - 0.5) * ddy))
            push!(hdiag, h)
        end
    end
    xf = min(cfg.x_fine, cfg.re_hi)
    pave!(0.0, xf, cfg.dx_fine, cfg.dy_fine)
    xf < cfg.re_hi && pave!(xf, cfg.re_hi, cfg.dx_coarse, cfg.dy_coarse)
    return _CSICover(centers, hdiag, fill(Inf, length(centers)))
end

_csi_covered(c::_CSICover) = isempty(c.centers)

# Register a shift and its disk(s): drop now-covered cells, update running shift distances.
# Radii of 0 are legal (shift contributes distance information but no coverage).
function _csi_cover!(c::_CSICover, σ::ComplexF64, r::Float64, r2::Float64 = 0.0)
    rmax = max(r, r2)
    keep = 0
    @inbounds for j in eachindex(c.centers)
        d = abs(c.centers[j] - σ)
        d + c.hdiag[j] ≤ rmax && continue
        keep += 1
        c.centers[keep] = c.centers[j]
        c.hdiag[keep] = c.hdiag[j]
        c.dmin[keep] = min(c.dmin[j], d)
    end
    resize!(c.centers, keep); resize!(c.hdiag, keep); resize!(c.dmin, keep)
    return c
end

# Next shift target: the uncovered cell farthest from all previous shifts (max-min distance
# spreads shifts instead of clumping at one hole). Returns nothing when covered.
function _csi_next_target(c::_CSICover, prev::Vector{ComplexF64})
    isempty(c.centers) && return nothing
    j = argmax(c.dmin)
    best = c.centers[j]
    # A repeat of an earlier (failed) shift would repeat the failure; nudge deterministically.
    if c.dmin[j] < 1e-3
        best += complex(0.011, 0.017) * (1 + length(prev) % 3)
    end
    return best
end

# Process one shift's results for one pencil: collect enumerated Ritz pairs (λ, residual),
# then extend coverage by the enumerated-trust disk.
function _csi_absorb!(cands::Vector{Tuple{ComplexF64,Float64}}, cover::_CSICover,
                      σ::ComplexF64, λcol, rescol, cfg::CertifiedSIConfig)
    enum = Tuple{Float64,ComplexF64}[]         # (|λ-σ|, λ)
    dbad = Inf                                 # nearest non-enumerated Ritz distance (cap)
    dmin = Inf                                 # nearest Ritz distance of any quality
    for i in eachindex(λcol)
        λ = λcol[i]
        isfinite(real(λ)) || continue
        dmin = min(dmin, abs(λ - σ))
        if rescol[i] ≤ cfg.enum_tol
            push!(enum, (abs(λ - σ), λ))
            push!(cands, (λ, rescol[i]))
        else
            dbad = min(dbad, abs(λ - σ))
        end
    end
    # Emptiness disk: no true mode lies within empty_frac × (nearest Ritz distance).
    r_empty = isfinite(dmin) ? cfg.empty_frac * dmin : 0.0
    r_trust = 0.0
    if length(enum) ≥ cfg.min_cert
        sort!(enum; by = first)
        q = min(cfg.trust, length(enum))
        r_trust = min(enum[q][1], dbad)
    end
    _csi_cover!(cover, σ, r_empty, r_trust)
    return
end

# Inverse/Rayleigh iteration from an enumerated (but not yet certified) Ritz value; no starting
# vector is needed because the shift already sits within the mode's basin. Returns (λ, resid).
function _csi_refine(A::AbstractMatrix{ComplexF64}, B::AbstractMatrix{ComplexF64},
                     λ::ComplexF64, iters::Int, tol::Float64)
    n = size(A, 1)
    x = fill(ComplexF64(1 / sqrt(n)), n)
    r = Inf
    for it in 1:iters
        F = lu(A - λ * B; check = false)
        issuccess(F) || return λ, r            # exact hit: λ is (numerically) an eigenvalue
        x = F \ (B * x)
        x ./= norm(x)
        Ax = A * x; Bx = B * x
        it > 1 && (λ = (x' * Ax) / (x' * Bx))  # hold λ once: lock onto the nearest mode first
        r = norm(Ax .- λ .* Bx) / max(norm(Ax) + abs(λ) * norm(Bx), eps())
        r ≤ tol && break
    end
    return λ, r
end

# Finalize one pencil: keep certified values as-is; refine relevant uncertified candidates
# (real(λ) > cert_lo) to full tolerance, flagging the pencil when any resists certification.
function _csi_finalize(cands::Vector{Tuple{ComplexF64,Float64}}, A, B,
                       cfg::CertifiedSIConfig)
    cert = ComplexF64[]
    todo = Tuple{ComplexF64,Float64}[]
    for (λ, res) in sort(cands; by = last)
        if res ≤ cfg.resid_tol
            push!(cert, λ)
        elseif real(λ) > cfg.cert_lo
            # pre-dedup refine inputs at the enumeration accuracy scale (cross-shift duplicates)
            any(t -> abs(t[1] - λ) < 1e-4, todo) && continue
            any(c -> abs(c - λ) < 1e-4, cert) && continue
            push!(todo, (λ, res))
        end
    end
    # Too many uncertified relevant candidates: one dense geev is cheaper and exact than
    # dozens of n³ LU refinements, so flag instead (bounds host cost on wide/crowded windows).
    length(todo) > cfg.max_refine && return _dedup_ritz(cert, cfg.dedup_tol), true
    flag = false
    for (λ0, _) in todo
        λr, rr = _csi_refine(A, B, λ0, cfg.refine, cfg.resid_tol)
        if rr ≤ cfg.resid_tol && abs(λr - λ0) ≤ cfg.max_move
            push!(cert, λr)
        else
            flag = true
        end
    end
    return _dedup_ritz(cert, cfg.dedup_tol), flag
end

# Base shifts: near-axis row dense inside the crowd band |ω| ≤ row_dense_w, row2_dy-spaced
# wings beyond it, plus the sparse mid-window row (full width).
function _csi_base_shifts(cfg::CertifiedSIConfig)
    w1 = min(cfg.im_max, cfg.row_dense_w)
    n1 = max(2, round(Int, 2 * w1 / cfg.row_dy) + 1)
    base = [complex(cfg.row_re, y) for y in range(-w1, w1; length = n1)]
    if cfg.im_max > w1 + 1e-9
        nw = ceil(Int, (cfg.im_max - w1) / cfg.row2_dy)
        for k in 1:nw, s in (-1.0, 1.0)
            push!(base, complex(cfg.row_re, s * (w1 + k * (cfg.im_max - w1) / nw)))
        end
    end
    nrow2 = max(2, round(Int, 2 * cfg.im_max / cfg.row2_dy) + 1)
    append!(base, complex(cfg.row2_re, y)
            for y in range(-cfg.im_max, cfg.im_max; length = nrow2))
    return base
end

"""
    certified_si_eigvals_batch(As, Bs; cfg=CertifiedSIConfig(), use_gpu=false,
                               dense_fallback=:cpu)
        -> (vals::Vector{Vector{ComplexF64}}, flagged::BitVector, nshifts::Vector{Int},
            reasons::Vector{Symbol})

Residual-certified eigenvalues of each pencil `(As[p], Bs[p])` near/inside the unstable window,
with a per-pencil geometric COVERAGE certificate: the union of certified shift-invert disks must
cover `[0,re_hi] × [-im_max,im_max]`; uncovered pencils receive adaptive per-pencil shifts
(batched rounds, compacted onto the still-uncovered subset) and, if still uncovered after
`max_rounds`, are `flagged` and — unless `dense_fallback === :none` — replaced by the exact
dense spectrum. `nshifts[p]` reports the shifts spent on pencil p; `reasons[p]` is `:ok`,
`:uncovered` (coverage budget exhausted), `:uncert` (a relevant candidate resisted refinement),
or `:both`.
"""
function certified_si_eigvals_batch(As::AbstractVector{<:AbstractMatrix{ComplexF64}},
                                    Bs::AbstractVector{<:AbstractMatrix{ComplexF64}};
                                    cfg::CertifiedSIConfig = CertifiedSIConfig(),
                                    use_gpu::Bool = false, dense_fallback::Symbol = :cpu)
    length(As) == length(Bs) || throw(DimensionMismatch("As and Bs must have equal length"))
    P = length(As)
    P == 0 && return Vector{Vector{ComplexF64}}(), falses(0), Int[], Symbol[]
    n = size(As[1], 1)
    all(A -> size(A, 1) == n, As) || throw(DimensionMismatch("all pencils must share one size n"))

    vals = Vector{Vector{ComplexF64}}(undef, P)
    flagged = falses(P)
    nshifts = zeros(Int, P)
    reasons = fill(:ok, P)

    # Shard so the GPU stack (A3+B3+G + n×M blocks) fits; CPU path takes one shard.
    maxP = (use_gpu && _CUDA_CSI_MAXP[] !== nothing) ? max(1, _CUDA_CSI_MAXP[](n, cfg.M)) : P
    for lo in 1:maxP:P
        rng = lo:min(lo + maxP - 1, P)
        _csi_solve_shard!(vals, flagged, reasons, nshifts, view(As, rng), view(Bs, rng), rng,
                          cfg, use_gpu)
    end

    _dense_fallback!(vals, flagged, As, Bs, dense_fallback)
    return vals, flagged, nshifts, reasons
end

function _csi_solve_shard!(vals, flagged, reasons, nshifts, As, Bs, out_rng,
                           cfg::CertifiedSIConfig, use_gpu::Bool)
    P = length(As)
    opener = (use_gpu && _CUDA_CSI_OPEN[] !== nothing) ? _CUDA_CSI_OPEN[] : _csi_open_cpu
    cands = [Tuple{ComplexF64,Float64}[] for _ in 1:P]
    covers = [_CSICover(cfg) for _ in 1:P]
    used = [ComplexF64[] for _ in 1:P]

    h = opener(As, Bs, cfg.M, cfg.Q)
    hidx = collect(1:P)                       # shard-local pencil index behind each handle lane
    try
        for σ in _csi_base_shifts(cfg)
            λ, res = h.solve(fill(σ, P))
            Threads.@threads for p in 1:P
                push!(used[p], σ)
                _csi_absorb!(cands[p], covers[p], σ, view(λ, :, p), view(res, :, p), cfg)
            end
        end

        active = collect(1:P)
        for _ in 1:cfg.max_rounds
            tg = Vector{Union{Nothing,ComplexF64}}(undef, length(active))
            Threads.@threads for i in eachindex(active)
                tg[i] = _csi_next_target(covers[active[i]], used[active[i]])
            end
            keep = findall(t -> t !== nothing, tg)
            isempty(keep) && break
            targets = Dict{Int,ComplexF64}(active[i] => ComplexF64(tg[i]) for i in keep)
            active = active[keep]
            # Compact the batched handle onto the still-uncovered pencils once wasted lanes
            # dominate (a reopen re-stacks the subset; it pays for itself within a round or two).
            if length(active) < length(hidx) ÷ 2
                h.close()
                h = opener(view(As, active), view(Bs, active), cfg.M, cfg.Q)
                hidx = copy(active)
            end
            σh = [haskey(targets, hidx[j]) ? targets[hidx[j]] : used[hidx[j]][end]
                  for j in eachindex(hidx)]
            λ, res = h.solve(σh)
            Threads.@threads for j in eachindex(hidx)
                p = hidx[j]
                haskey(targets, p) || continue
                push!(used[p], targets[p])
                _csi_absorb!(cands[p], covers[p], targets[p], view(λ, :, j), view(res, :, j),
                             cfg)
            end
        end

        Threads.@threads for p in 1:P
            po = out_rng[p]
            uncovered = !_csi_covered(covers[p])
            nshifts[po] = length(used[p])
            vals[po], uncert = _csi_finalize(cands[p], As[p], Bs[p], cfg)
            flagged[po] = uncovered | uncert
            reasons[po] = uncovered ? (uncert ? :both : :uncovered) :
                          (uncert ? :uncert : :ok)
        end
    finally
        h.close()
    end
    return
end

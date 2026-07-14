using Test
using TJLF
using LinearAlgebra
using Random

# Batched shift-invert dominant-mode eigensolver (TJLF.si_eigvals_batch). Tests the CPU reference
# path (GPU path is exercised separately on a CUDA node). We verify that every dense eigenvalue
# lying near a configured shift is recovered by the shift-invert solver.
@testset "batched shift-invert eigensolver" begin
    Random.seed!(42)
    n = 60
    P = 5
    shifts = ComplexF64[0.0, 0.5im, -0.5im, 0.3 + 0.0im]
    cfg = TJLF.SIConfig(shifts = shifts, M = 12, Q = 20)

    As = Vector{Matrix{ComplexF64}}(undef, P)
    Bs = Vector{Matrix{ComplexF64}}(undef, P)
    targets = Vector{Vector{ComplexF64}}(undef, P)
    for p in 1:P
        # Controlled, well-separated spectrum: two isolated eigenvalues near each shift plus a
        # distant damped bulk (|lambda| >= 3). Normal matrix A = Q diag(lambda) Q^H keeps the
        # eigenvalues exactly `lambda` and perfectly conditioned, so shift-invert must recover
        # the near-shift eigenvalues tightly. B = I -> standard problem M = A.
        near = ComplexF64[]
        for σ in shifts
            push!(near, σ + 0.03 * (randn() + im * randn()))
            push!(near, σ + 0.06 * (randn() + im * randn()))
        end
        bulk = ComplexF64[(3 + 2 * rand()) * cis(2π * rand()) for _ in 1:(n - length(near))]
        λ = vcat(near, bulk)
        Q = Matrix(qr(randn(ComplexF64, n, n)).Q)
        As[p] = Q * Diagonal(λ) * Q'
        Bs[p] = Matrix{ComplexF64}(I, n, n)
        targets[p] = near
    end

    got = TJLF.si_eigvals_batch(As, Bs; cfg = cfg, use_gpu = false)
    @test length(got) == P

    for p in 1:P
        # Each planted near-shift eigenvalue must be recovered accurately.
        for λ in targets[p]
            err = minimum(abs(λ - g) for g in got[p]; init = Inf)
            @test err < 1e-5
        end
        @test !isempty(got[p])
    end
end

# Contour-integral (Sakurai-Sugiura Rayleigh-Ritz) batched eigensolver
# (TJLF.contour_eigvals_batch). CPU moment path; the GPU path computes the identical moments
# with batched cuBLAS and is exercised on a CUDA node. The contract is stronger than
# shift-invert's: EVERY eigenvalue inside the rectangular window must be recovered
# (residual-certified) and everything outside excluded.
@testset "contour-integral batched eigensolver" begin
    Random.seed!(7)
    n = 80
    P = 4
    cfg = TJLF.ContourConfig(re_lo = -0.05, re_hi = 0.4, im_max = 1.2,
                             n_long = 12, n_short = 6, L = 16, K = 2)

    As = Vector{Matrix{ComplexF64}}(undef, P)
    Bs = Vector{Matrix{ComplexF64}}(undef, P)
    inside = Vector{Vector{ComplexF64}}(undef, P)
    for p in 1:P
        # Planted spectrum: eigenvalues scattered THROUGHOUT the window including the corners
        # (small gamma, large |freq| — exactly where fixed-shift SI failed), a few weakly damped
        # modes just outside the window, and a strongly damped bulk. Generalized pencil with a
        # well-conditioned random B.
        λin = [complex(-0.03 + 0.4 * rand(), -1.1 + 2.2 * rand()) for _ in 1:(5 + p)]
        push!(λin, complex(0.01, 1.1), complex(0.35, -1.1))           # corner modes
        λin = [λ for λ in λin if TJLF._contour_inside(λ, cfg) &&
               real(λ) > cfg.re_lo + 0.01 && abs(imag(λ)) < cfg.im_max - 0.05]
        λout = ComplexF64[-0.3 + 0.3im, -0.3 - 0.9im, 0.1 + 1.6im]    # outside, near-ish window
        bulk = ComplexF64[(20 + 30 * rand()) * cis(π * (0.55 + 0.9 * rand())) for
                          _ in 1:(n - length(λin) - length(λout))]
        λ = vcat(λin, λout, bulk)
        Q = Matrix(qr(randn(ComplexF64, n, n)).Q)
        M = Q * Diagonal(λ) * Q'
        Bs[p] = Matrix{ComplexF64}(I, n, n) + 0.1 * randn(ComplexF64, n, n)
        As[p] = Bs[p] * M                       # B⁻¹A = M has spectrum exactly λ
        inside[p] = λin
    end

    vals, flagged, ranks = TJLF.contour_eigvals_batch(As, Bs; cfg = cfg, use_gpu = false,
                                                      dense_fallback = :none)
    @test length(vals) == P
    for p in 1:P
        @test !flagged[p]
        # Every interior eigenvalue recovered tightly (the anti-"coverage gap" contract).
        for λ in inside[p]
            @test minimum(abs(λ - g) for g in vals[p]; init = Inf) < 1e-8
        end
        # No spurious eigenvalues: everything returned matches a planted interior mode.
        for g in vals[p]
            @test minimum(abs(λ - g) for λ in inside[p]; init = Inf) < 1e-8
        end
        @test ranks[p] ≥ length(inside[p])
    end

    # Capacity saturation must be detected and (with dense_fallback=:cpu) repaired exactly:
    # plant more interior eigenvalues than the 2K·L subspace capacity.
    cfg_small = TJLF.ContourConfig(re_lo = -0.05, re_hi = 0.4, im_max = 1.2,
                                   n_long = 12, n_short = 6, L = 3, K = 1)
    λin = [complex(-0.02 + 0.35 * rand(), -1.0 + 2.0 * rand()) for _ in 1:8]
    bulk = ComplexF64[(20 + 30 * rand()) * cis(π * (0.55 + 0.9 * rand())) for _ in 1:(n - 8)]
    Q = Matrix(qr(randn(ComplexF64, n, n)).Q)
    A1 = Q * Diagonal(vcat(λin, bulk)) * Q'
    B1 = Matrix{ComplexF64}(I, n, n)
    _, fl, _ = TJLF.contour_eigvals_batch([A1], [B1]; cfg = cfg_small, use_gpu = false,
                                          dense_fallback = :none)
    @test fl[1]
    v2, fl2, _ = TJLF.contour_eigvals_batch([A1], [B1]; cfg = cfg_small, use_gpu = false,
                                            dense_fallback = :cpu)
    @test fl2[1]
    @test length(v2[1]) == n          # dense fallback returns the full spectrum
    for λ in λin
        @test minimum(abs(λ - g) for g in v2[1]; init = Inf) < 1e-8
    end
end

# Coverage-certified adaptive multi-shift SI (TJLF.certified_si_eigvals_batch). CPU solve path;
# the GPU handle computes the same quantities with batched cuBLAS and is exercised on a CUDA
# node. Planted spectra mimic the REAL TJLF structure that killed both the fixed-shift SI
# (coverage gaps) and the contour solver (no spectral gap): a dense crowd of weakly damped
# modes hugging gamma <= 0, unstable leaders scattered across the window including far corners,
# and a strongly damped bulk.
@testset "coverage-certified adaptive SI" begin
    Random.seed!(21)
    n = 120
    P = 4
    cfg = TJLF.CertifiedSIConfig(re_hi = 0.6, im_max = 1.5, row_dy = 0.5,
                                 M = 16, Q = 16, trust = 8, max_rounds = 60)

    As = Vector{Matrix{ComplexF64}}(undef, P)
    Bs = Vector{Matrix{ComplexF64}}(undef, P)
    unstable = Vector{Vector{ComplexF64}}(undef, P)
    for p in 1:P
        # Unstable modes anywhere in the window — including spots far from the initial near-axis
        # shift row (high-drive corners), the historical fixed-shift blind spot.
        λu = [complex(0.02 + 0.55 * rand(), -1.4 + 2.8 * rand()) for _ in 1:(2 + p)]
        crowd = [complex(-0.06 * rand(), -1.4 + 2.8 * rand()) for _ in 1:40]
        bulk = ComplexF64[(15 + 25 * rand()) * cis(π * (0.55 + 0.9 * rand()))
                          for _ in 1:(n - length(λu) - length(crowd))]
        Q = Matrix(qr(randn(ComplexF64, n, n)).Q)
        M = Q * Diagonal(vcat(λu, crowd, bulk)) * Q'
        Bs[p] = Matrix{ComplexF64}(I, n, n) + 0.1 * randn(ComplexF64, n, n)
        As[p] = Bs[p] * M
        unstable[p] = λu
    end

    vals, flagged, nshifts, reasons = TJLF.certified_si_eigvals_batch(As, Bs; cfg = cfg,
                                                                      use_gpu = false,
                                                                      dense_fallback = :none)
    @test all(r -> r in (:ok, :uncovered, :uncert, :both), reasons)
    @test length(vals) == P
    for p in 1:P
        if !flagged[p]
            # Coverage certified -> EVERY planted unstable mode must be present and tight.
            for λ in unstable[p]
                @test minimum(abs(λ - g) for g in vals[p]; init = Inf) < 1e-7
            end
        end
        # Everything returned is certified -> must match a true eigenvalue of the pencil.
        Λtrue = eigvals(Bs[p] \ As[p])
        for g in vals[p]
            @test minimum(abs(λ - g) for λ in Λtrue; init = Inf) < 1e-6
        end
        @test nshifts[p] ≥ 2
    end
    # The window must be coverable for these well-separated planted spectra.
    @test count(flagged) == 0

    # Starved config (tiny round budget) must flag, and dense fallback must repair exactly.
    cfg_starved = TJLF.CertifiedSIConfig(re_hi = 0.6, im_max = 1.5, row_dy = 3.0,
                                         M = 6, Q = 4, trust = 3, min_cert = 2, max_rounds = 0)
    v3, fl3, _ = TJLF.certified_si_eigvals_batch(As[1:1], Bs[1:1]; cfg = cfg_starved,
                                                 use_gpu = false, dense_fallback = :cpu)
    @test fl3[1]
    @test length(v3[1]) == n
    for λ in unstable[1]
        @test minimum(abs(λ - g) for g in v3[1]; init = Inf) < 1e-8
    end
end

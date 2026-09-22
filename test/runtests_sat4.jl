using Test
using TJLF
using ForwardDiff  # provided by test/Project.toml; triggers TJLFForwardDiffExt

# SAT4 = SAT0 intensity formula on the SAT2/3 linear physics, with calibratable coefficients.
# Deck: tglf01 (GA standard, NS=2, VEXB_SHEAR=0 -> single pass) with SAT_RULE overridden to 4.

const SAT4_DECK = joinpath(@__DIR__, "tglf_regression", "tglf01", "input.tglf")
const SAT2_DECK = joinpath(@__DIR__, "test_SAT_rules", "SAT2", "input.tglf")

function sat4_input(deck=SAT4_DECK; sat_rule=4)
    inp = readInput(deck)
    inp.SAT_RULE = sat_rule
    TJLF.apply_presets!(inp)
    return inp
end

# cached linear solve (same split as runtests_core.jl); the saturation coefficients only
# enter sum_ky_spectrum, so the cache can be re-used for every coefficient vector
function sat4_linear_cache(inp)
    satParams = get_sat_params(inp)
    inp.KY_SPECTRUM .= get_ky_spectrum(inp, satParams.grad_r0)
    hermite = gauss_hermite(inp)
    tm = tjlf_TM(inp, satParams, hermite)
    single_pass = inp.ALPHA_QUENCH != 0 || inp.VEXB_SHEAR * inp.SIGN_IT == 0.0
    gamma = single_pass ? tm.firstPass_eigenvalue[:, :, 1] : tm.secondPass_eigenvalue[:, :, 1]
    vzf, kymax, jmax = TJLF.get_zonal_mixing(inp, satParams, tm.firstPass_eigenvalue[1, :, 1])
    return (; satParams, gamma, QL=tm.QL_weights, vzf, kymax, jmax)
end

sat4_fluxes(inp, c) = sum_ky_spectrum(inp, c.satParams, c.gamma, c.QL;
                                      vzf_out_param=c.vzf, kymax_out_param=c.kymax, jmax_out_param=c.jmax)[1]

# independent transcription of tglf_LS.f90 get_intensity (igeo=1) with free constants
function sat4_reference(inp, sp, gamma, kx0_e; C_NORM=inp.C_NORM, C_EXP=inp.C_EXP,
                        C_COEFF=inp.C_COEFF, C_ETG=inp.C_ETG)
    ky = inp.KY_SPECTRUM
    out = zeros(length(ky), inp.NMODES)
    pol = sum(inp.ZS[s]^2 * inp.AS[s] / inp.TAUS[s] for s in 1:inp.NS)
    pols = (pol / abs(inp.AS[1] * inp.ZS[1]^2))^2
    measure = sqrt(inp.TAUS[1] * inp.MASS[2])
    for j in eachindex(ky), i in 1:inp.NMODES
        g = gamma[i, j]
        g > 0 || continue
        ks = ky[j] * measure / abs(inp.ZS[1])
        cnorm = C_NORM * pols
        ks > 1 && (cnorm /= ks^C_ETG)
        wd0 = ks * sqrt(inp.TAUS[1] / inp.MASS[2]) / sp.R_unit
        gnet = g / wd0
        I = cnorm * wd0^2 * (gnet^C_EXP + C_COEFF * gnet) / ky[j]^4
        if inp.ALPHA_QUENCH == 0 && abs(kx0_e[j]) > 0
            I /= (1 + 0.56 * kx0_e[j]^2)^2
            I /= (1 + (1.15 * kx0_e[j])^4)^2
        end
        out[j, i] = I * sp.SAT_geo0 * measure
    end
    return out
end

@testset "SAT4 saturation rule" begin

    @testset "presets and input validation" begin
        inp = sat4_input()
        @test inp.SAT_RULE == 4
        @test inp.UNITS == "CGYRO"
        @test inp.XNU_MODEL == 3
        @test inp.WDIA_TRAPPED == 1.0
        @test isnan(inp.C_B)               # NaN sentinel = Fortran constant, allowed by checkInput
        @test inp.SIG_B == 0.34
        @test inp.BOUNCE_COEFF == 3.0
        @test (inp.C_NORM, inp.C_EXP, inp.C_COEFF, inp.C_ETG) == (1.82770384, 1.39786897, 0.36017009, 1.25)
        TJLF.checkInput(inp)               # must not throw
        bad = readInput(SAT4_DECK); bad.SAT_RULE = 5
        @test_throws AssertionError TJLF.checkInput(bad)
    end

    @testset "input.tglf round trip of the TJLF-only keys" begin
        inp = sat4_input()
        inp.C_NORM = 2.5; inp.C_EXP = 1.2; inp.C_COEFF = 0.1; inp.C_ETG = 1.4
        inp.C_B = 0.4; inp.SIG_B = 0.3; inp.BOUNCE_COEFF = 2.0
        path = tempname() * ".tglf"
        TJLF.save(inp, path)
        back = readInput(path)
        for k in (:C_NORM, :C_EXP, :C_COEFF, :C_ETG, :C_B, :SIG_B, :BOUNCE_COEFF)
            @test getfield(back, k) == getfield(inp, k)
        end
        @test back.SAT_RULE == 4 && back.UNITS == "CGYRO" && back.XNU_MODEL == 3
        rm(path)
    end

    @testset "intensity matches the SAT0 formula (tglf01)" begin
        inp = sat4_input()
        c = sat4_linear_cache(inp)
        res = TJLF.intensity_sat(inp, c.satParams, c.gamma, c.QL, 2.0, true;
                                 vzf_out_param=c.vzf, kymax_out_param=c.kymax, jmax_out_param=c.jmax)
        ref = sat4_reference(inp, c.satParams, c.gamma, res.kx0_e)
        @test size(res.phinorm) == (length(inp.KY_SPECTRUM), inp.NMODES)
        @test all(isfinite, res.phinorm)
        @test any(>(0), res.phinorm)
        @test isapprox(res.phinorm, ref; rtol=1e-12, atol=0.0)
        # SAT4 has no QLA correction
        phinorm, QLA_P, QLA_E, QLA_O = TJLF.intensity_sat(inp, c.satParams, c.gamma, c.QL;
                                 vzf_out_param=c.vzf, kymax_out_param=c.kymax, jmax_out_param=c.jmax)
        @test phinorm == res.phinorm
        @test all(==(1.0), QLA_P) && all(==(1.0), QLA_E) && all(==(1.0), QLA_O)
    end

    @testset "coefficient scaling laws" begin
        inp = sat4_input()
        c = sat4_linear_cache(inp)
        f0 = sat4_fluxes(inp, c)
        @test size(f0) == (3, inp.NS, 5)
        @test TJLF.Qe(f0) > 0
        inp.C_NORM *= 2
        f2 = sat4_fluxes(inp, c)
        @test isapprox(f2, 2 .* f0; rtol=1e-13)   # intensity is linear in C_NORM
        inp.C_NORM /= 2
        # C_EXP=1, C_COEFF=0 -> intensity proportional to gamma at fixed ky
        inp.C_EXP = 1.0; inp.C_COEFF = 0.0
        res = TJLF.intensity_sat(inp, c.satParams, c.gamma, c.QL, 2.0, true;
                                 vzf_out_param=c.vzf, kymax_out_param=c.kymax, jmax_out_param=c.jmax)
        for j in eachindex(inp.KY_SPECTRUM)
            unstable = [i for i in 1:inp.NMODES if c.gamma[i, j] > 0]
            length(unstable) >= 2 || continue
            r = [res.phinorm[j, i] / c.gamma[i, j] for i in unstable]
            @test all(isapprox.(r, r[1]; rtol=1e-12))
        end
    end

    @testset "golden regression: run_tjlf SAT4 on tglf01" begin
        # TJLF self-consistency baseline captured on the sat4 branch (master 836e31a + SAT4 port);
        # layout of sum(flux; dims=1) flattened = [G_e, G_i, Q_e, Q_i, Πtor_e, Πtor_i, Πpar_e, Πpar_i, X_e, X_i]
        inp = sat4_input()
        fl = vec(sum(TJLF.run_tjlf(inp); dims=1))
        golden = [-1.1384094678645216, -1.138409467864521, 13.465180684019883, 30.778337829532077,
                  0.0, 0.0, 0.0, 0.0, 3.955288680090136, -3.955288680090139]
        for i in eachindex(golden)
            @test isapprox(fl[i], golden[i]; rtol=5e-3, atol=1e-8)
        end
        # ordering of the rules on this deck: SAT4 default sits below SAT0 (Qe 16.2) and SAT2 (22.8)
        @test 10 < fl[3] < 16
    end

    @testset "defaults reproduce the unmodified collision/trapping model (SAT2 deck)" begin
        a = readInput(SAT2_DECK)
        fa = TJLF.run_tjlf(a)
        b = readInput(SAT2_DECK)
        b.C_B = 0.315; b.SIG_B = 0.34; b.BOUNCE_COEFF = 3.0   # the Fortran literals, explicitly
        fb = TJLF.run_tjlf(b)
        @test fa == fb
        # SAT4 also runs on an NS=3, electromagnetic, sheared deck
        d = sat4_input(SAT2_DECK)
        fd = TJLF.run_tjlf(d)
        @test all(isfinite, fd)
        @test TJLF.Qe(fd) > 0
    end

    @testset "ForwardDiff through sum_ky_spectrum" begin
        @test !isnothing(Base.get_extension(TJLF, :TJLFForwardDiffExt))
        inp = sat4_input()
        c = sat4_linear_cache(inp)

        function promote_input(::Type{T}, base::TJLF.InputTJLF{Float64}) where {T<:Real}
            out = TJLF.InputTJLF{T}(base.NS, length(base.KY_SPECTRUM))
            for fn in fieldnames(TJLF.InputTJLF)
                v = getfield(base, fn)
                if v isa Float64
                    setfield!(out, fn, T(v))
                elseif v isa Vector{Float64}
                    setfield!(out, fn, T.(v))
                elseif v isa Vector{ComplexF64}
                    setfield!(out, fn, Complex{T}.(v))
                else
                    setfield!(out, fn, v)
                end
            end
            return out
        end
        promote_val(::Type{T}, v::Float64) where {T} = T(v)
        promote_val(::Type{T}, v::AbstractArray{Float64}) where {T} = T.(v)
        promote_val(::Type{T}, v) where {T} = v
        function promote_satparams(::Type{T}, sp::TJLF.SaturationParameters{Float64}) where {T}
            vals = map(fn -> promote_val(T, getfield(sp, fn)), fieldnames(TJLF.SaturationParameters))
            return TJLF.SaturationParameters{T}(vals...)
        end

        function Qe_of(θ::AbstractVector{T}) where {T<:Real}
            i = promote_input(T, inp)
            i.C_NORM, i.C_EXP, i.C_COEFF = θ[1], θ[2], θ[3]
            sp = promote_satparams(T, c.satParams)
            fl = sum_ky_spectrum(i, sp, T.(c.gamma), T.(c.QL);
                                 vzf_out_param=T(c.vzf), kymax_out_param=T(c.kymax), jmax_out_param=c.jmax)[1]
            return TJLF.Qe(fl)
        end

        θ0 = [inp.C_NORM, inp.C_EXP, inp.C_COEFF]
        q0 = Qe_of(θ0)
        g = ForwardDiff.gradient(Qe_of, θ0)
        @test all(isfinite, g)
        @test isapprox(g[1], q0 / θ0[1]; rtol=1e-10)        # exact by linearity in C_NORM
        h = 1e-6
        for k in 2:3
            θp = copy(θ0); θp[k] += h
            θm = copy(θ0); θm[k] -= h
            fd = (Qe_of(θp) - Qe_of(θm)) / (2h)
            @test isapprox(g[k], fd; rtol=1e-5)
        end
    end
end

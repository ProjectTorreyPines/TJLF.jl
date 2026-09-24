using Test
using TJLF
using ForwardDiff  # provided by test/Project.toml; triggers TJLFForwardDiffExt

# SAT1/SAT2/SAT3 saturation-rule constants exposed as InputTJLF fields (SAT1_*, SAT2_*, SAT3_*).
# Their defaults must reproduce the previously hard-coded Fortran literals bit-exactly.

const SATC_DECKS = Dict(
    "SAT1" => joinpath(@__DIR__, "test_SAT_rules", "SAT1", "input.tglf"),
    "SAT2" => joinpath(@__DIR__, "test_SAT_rules", "SAT2", "input.tglf"),
    "SAT3" => joinpath(@__DIR__, "test_SAT_rules", "SAT3", "input.tglf"),
)
const SATC_TGLF01 = joinpath(@__DIR__, "tglf_regression", "tglf01", "input.tglf")

# independent transcription of the Fortran literals (tglf_multiscale_spectrum.f90 / intensity_sat)
const SATC_LITERALS = Dict(
    :SAT1_CNORM => 14.29, :SAT1_CZ1 => 0.48, :SAT1_CZ2 => 1.0, :SAT1_ETG_STREAMER => 1.05, :SAT1_CKY => 3.0,
    :SAT1_AX => 1.15, :SAT1_AY => 0.56,
    :SAT2_B0 => 0.76, :SAT2_B1 => 1.22, :SAT2_B2 => 3.55, :SAT2_B2_SINGLE => 3.74, :SAT2_B3 => 1.0,
    :SAT2_CZ2 => 1.05, :SAT2_CKY => 3.0, :SAT2_AX => 1.21, :SAT2_AY => 1.0, :SAT2_KYETG => 1000.0,
    :SAT3_Y_ITG => 3.3, :SAT3_Y_TEM => 12.7, :SAT3_SCAL => 0.82, :SAT3_KMIN => 0.685, :SAT3_COVERB => -0.751,
    :SAT3_C1 => -2.42, :SAT3_K0 => 0.6, :SAT3_KP => 2.0, :SAT3_X_ITG => 0.8, :SAT3_X_TEM => 1.0,
    :SAT3_QLA_P_ITG => 1.1, :SAT3_QLA_P_TEM => 0.6, :SAT3_QLA_E_ITG => 0.75, :SAT3_QLA_E_TEM => 0.6, :SAT3_QLA_O => 0.8,
)

# (Qe, Qi, Γe) captured with the hard-coded literals (sat4 branch 56c91b9) on this machine
const SATC_GOLDEN = Dict(
    "SAT1" => (2.5980298355708413, 5.097402979045156, 0.22393062392652113),
    "SAT2" => (5.352639266058449, 15.291135349661264, 0.6002988303946646),
    "SAT3" => (3.873989360845364, 9.809559273004322, -0.04184899128611297),
    "tglf01_sat1" => (16.857891666690318, 44.02270123034034, -2.249128250763428),
    "tglf01_sat2" => (22.763940659694555, 65.27557048777821, -3.0235301405715145),
    "tglf01_sat3" => (26.75900850370548, 78.38902803825452, -5.490624520158389),
    "SAT1_cgyro" => (5.7970688384686495, 13.709119666179783, 0.6999910892875976),
)

function satc_input(deck; sat_rule=nothing, units=nothing)
    inp = readInput(deck)
    sat_rule === nothing || (inp.SAT_RULE = sat_rule)
    units === nothing || (inp.UNITS = units)
    TJLF.apply_presets!(inp)
    return inp
end

satc_triple(fl) = (TJLF.Qe(fl), TJLF.Qi(fl), TJLF.Γe(fl))

function satc_linear_cache(inp)
    satParams = get_sat_params(inp)
    inp.KY_SPECTRUM .= get_ky_spectrum(inp, satParams.grad_r0)
    hermite = gauss_hermite(inp)
    tm = tjlf_TM(inp, satParams, hermite)
    single_pass = inp.ALPHA_QUENCH != 0 || inp.VEXB_SHEAR * inp.SIGN_IT == 0.0
    gamma = single_pass ? tm.firstPass_eigenvalue[:, :, 1] : tm.secondPass_eigenvalue[:, :, 1]
    vzf, kymax, jmax = TJLF.get_zonal_mixing(inp, satParams, tm.firstPass_eigenvalue[1, :, 1])
    return (; satParams, gamma, QL=tm.QL_weights, vzf, kymax, jmax)
end

function satc_promote_input(::Type{T}, base::TJLF.InputTJLF{Float64}) where {T<:Real}
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
satc_promote_val(::Type{T}, v::Float64) where {T} = T(v)
satc_promote_val(::Type{T}, v::AbstractArray{Float64}) where {T} = T.(v)
satc_promote_val(::Type{T}, v) where {T} = v
function satc_promote_satparams(::Type{T}, sp::TJLF.SaturationParameters{Float64}) where {T}
    vals = map(fn -> satc_promote_val(T, getfield(sp, fn)), fieldnames(TJLF.SaturationParameters))
    return TJLF.SaturationParameters{T}(vals...)
end

# fluxes from a cached linear solve; SAT1 (as run_tjlf) gets no zonal-mixing kwargs
function satc_fluxes(inp::TJLF.InputTJLF{T}, c) where {T}
    sp = satc_promote_satparams(T, c.satParams)
    if inp.SAT_RULE == 1
        return sum_ky_spectrum(inp, sp, T.(c.gamma), T.(c.QL))[1]
    end
    return sum_ky_spectrum(inp, sp, T.(c.gamma), T.(c.QL);
                           vzf_out_param=T(c.vzf), kymax_out_param=T(c.kymax), jmax_out_param=c.jmax)[1]
end

@testset "SAT1/2/3 constants as inputs" begin
    @testset "defaults equal the Fortran literals" begin
        inp = TJLF.InputTJLF{Float64}(2, 3)
        for (k, v) in SATC_LITERALS
            @test getfield(inp, k) == v
        end
        @test all(k -> hasfield(TJLF.InputTJLF, k), TJLF.TJLF_ONLY_KEYS)
        @test all(k -> k in TJLF.TJLF_ONLY_KEYS, keys(SATC_LITERALS))
        @test all(k -> isnan(getfield(TJLF.InputTGLF(), k)), keys(SATC_LITERALS))
    end

    @testset "defaults reproduce the hard-coded rules (bit-exact vs explicit literals, golden vs 56c91b9)" begin
        cases = [("SAT1", satc_input(SATC_DECKS["SAT1"])), ("SAT2", satc_input(SATC_DECKS["SAT2"])),
                 ("SAT3", satc_input(SATC_DECKS["SAT3"])), ("SAT1_cgyro", satc_input(SATC_DECKS["SAT1"]; units="CGYRO")),
                 ("tglf01_sat1", satc_input(SATC_TGLF01; sat_rule=1)), ("tglf01_sat2", satc_input(SATC_TGLF01; sat_rule=2)),
                 ("tglf01_sat3", satc_input(SATC_TGLF01; sat_rule=3))]
        for (name, a) in cases
            fa = TJLF.run_tjlf(a)
            b = deepcopy(a)
            for (k, v) in SATC_LITERALS
                setfield!(b, k, v)
            end
            @test fa == TJLF.run_tjlf(b)
            ta = satc_triple(fa)
            for i in 1:3
                @test isapprox(ta[i], SATC_GOLDEN[name][i]; rtol=1e-6, atol=1e-10)
            end
        end
    end

    @testset "save / readInput round trip" begin
        inp = satc_input(SATC_DECKS["SAT3"])
        for (i, k) in enumerate(sort(collect(keys(SATC_LITERALS))))
            setfield!(inp, k, getfield(inp, k) * (1 + 0.01 * i))
        end
        mktempdir() do dir
            path = joinpath(dir, "input.tglf")
            TJLF.save(inp, path)
            back = readInput(path)
            for k in keys(SATC_LITERALS)
                @test getfield(back, k) == getfield(inp, k)
            end
        end
    end

    @testset "scaling laws" begin
        # SAT1: flux ∝ SAT1_CNORM
        i1 = satc_input(SATC_DECKS["SAT1"]); c1 = satc_linear_cache(i1)
        f = satc_fluxes(i1, c1); i1.SAT1_CNORM *= 2
        @test isapprox(satc_fluxes(i1, c1), 2 .* f; rtol=1e-13)
        # SAT2: flux ∝ SAT2_B2 (NMODES > 1 here) and independent of SAT2_B2_SINGLE
        i2 = satc_input(SATC_DECKS["SAT2"]); c2 = satc_linear_cache(i2)
        @test i2.NMODES > 1
        f = satc_fluxes(i2, c2); i2.SAT2_B2_SINGLE = 99.0
        @test satc_fluxes(i2, c2) == f
        i2.SAT2_B2 *= 2
        @test isapprox(satc_fluxes(i2, c2), 2 .* f; rtol=1e-13)
        # SAT3: joint ×2 of (Y_ITG, Y_TEM, SCAL) doubles the flux (the connecting quadratic depends on Ys/YTs only)
        i3 = satc_input(SATC_DECKS["SAT3"]); c3 = satc_linear_cache(i3)
        f = satc_fluxes(i3, c3)
        i3.SAT3_Y_ITG *= 2; i3.SAT3_Y_TEM *= 2; i3.SAT3_SCAL *= 2
        @test isapprox(satc_fluxes(i3, c3), 2 .* f; rtol=1e-12)
        # SAT3 QLA_P scales the particle flux only
        i3 = satc_input(SATC_DECKS["SAT3"]); i3.KY_SPECTRUM .= get_ky_spectrum(i3, c3.satParams.grad_r0)
        i3.SAT3_QLA_P_ITG *= 3; i3.SAT3_QLA_P_TEM *= 3
        f3 = satc_fluxes(i3, c3)
        @test isapprox(f3[:, :, 1], 3 .* f[:, :, 1]; rtol=1e-13)
        @test f3[:, :, 2] == f[:, :, 2]
    end

    @testset "ForwardDiff through sum_ky_spectrum, per rule" begin
        @test !isnothing(Base.get_extension(TJLF, :TJLFForwardDiffExt))
        specs = [("SAT1", :SAT1_CNORM, :SAT1_CZ2), ("SAT2", :SAT2_B2, :SAT2_B0), ("SAT3", :SAT3_Y_ITG, :SAT3_KMIN)]
        for (name, level, shape) in specs
            inp = satc_input(SATC_DECKS[name])
            c = satc_linear_cache(inp)
            function Qe_of(θ::AbstractVector{T}) where {T<:Real}
                i = satc_promote_input(T, inp)
                setfield!(i, level, θ[1]); setfield!(i, shape, θ[2])
                return TJLF.Qe(satc_fluxes(i, c))
            end
            θ0 = [getfield(inp, level), getfield(inp, shape)]
            q0 = Qe_of(θ0)
            g = ForwardDiff.gradient(Qe_of, θ0)
            @test all(isfinite, g)
            if name != "SAT3"
                @test isapprox(g[1], q0 / θ0[1]; rtol=1e-10)   # exact by linearity in the level constant
            end
            h = 1e-5 * abs(θ0[2])
            fd = (Qe_of([θ0[1], θ0[2] + h]) - Qe_of([θ0[1], θ0[2] - h])) / (2h)
            @test isapprox(g[2], fd; rtol=1e-4, atol=1e-8 * abs(q0))
        end
    end
end

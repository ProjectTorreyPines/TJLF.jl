using Test
using TJLF
using ForwardDiff

@testset "ForwardDiff AD" begin

    # Helper: construct InputTJLF{T} from a Float64 baseline
    function convert_input_tjlf(::Type{T}, base::TJLF.InputTJLF{Float64}) where {T<:Real}
        ns = base.NS
        nky = length(base.KY_SPECTRUM)
        inp = TJLF.InputTJLF{T}(ns, nky)
        for fn in fieldnames(TJLF.InputTJLF)
            v = getfield(base, fn)
            ismissing(v) && continue
            if v isa Float64
                setfield!(inp, fn, T(v))
            elseif v isa Vector{Float64}
                setfield!(inp, fn, T.(v))
            elseif v isa Vector{ComplexF64}
                setfield!(inp, fn, Complex{T}.(v))
            else
                setfield!(inp, fn, v)
            end
        end
        return inp
    end

    # Use tglf01 with SAT_RULE=0 (XNUE=0.0, XNU_MODEL=2)
    # Exercises the cnuei=0 code path (guarded sqrt/pow in eigensolver)
    basedir = joinpath(@__DIR__, "tglf_regression", "tglf01")
    input_tjlf = TJLF.readInput(joinpath(basedir, "input.tglf"))

    # --- Test 1: ∂Qe/∂RLTS_2 ---
    @testset "∂Qe/∂RLTS_2" begin
        function tjlf_Qe(x::T) where {T<:Real}
            inp = convert_input_tjlf(T, input_tjlf)
            inp.RLTS[2] = x
            QL = TJLF.run_tjlf(inp)
            return TJLF.Qe(QL)
        end

        x0 = input_tjlf.RLTS[2]

        # Finite differences (reference)
        h = 1e-5
        dQe_fd = (tjlf_Qe(x0 + h) - tjlf_Qe(x0 - h)) / (2h)

        # ForwardDiff
        dQe_ad = ForwardDiff.derivative(tjlf_Qe, x0)

        @test isfinite(dQe_ad)
        @test isfinite(dQe_fd)
        relerr = abs(dQe_fd) > 0 ? abs(dQe_ad - dQe_fd) / abs(dQe_fd) : abs(dQe_ad - dQe_fd)
        @test relerr < 1e-4
    end

    # --- Test 2: ∂Qi/∂RLTS_2 ---
    @testset "∂Qi/∂RLTS_2" begin
        function tjlf_Qi(x::T) where {T<:Real}
            inp = convert_input_tjlf(T, input_tjlf)
            inp.RLTS[2] = x
            QL = TJLF.run_tjlf(inp)
            return TJLF.Qi(QL)
        end

        x0 = input_tjlf.RLTS[2]

        h = 1e-5
        dQi_fd = (tjlf_Qi(x0 + h) - tjlf_Qi(x0 - h)) / (2h)
        dQi_ad = ForwardDiff.derivative(tjlf_Qi, x0)

        @test isfinite(dQi_ad)
        @test isfinite(dQi_fd)
        relerr = abs(dQi_fd) > 0 ? abs(dQi_ad - dQi_fd) / abs(dQi_fd) : abs(dQi_ad - dQi_fd)
        @test relerr < 1e-4
    end

    # --- Test 3: ∂Qe/∂RLNS_2 (density gradient) ---
    @testset "∂Qe/∂RLNS_2" begin
        function tjlf_Qe_rlns(x::T) where {T<:Real}
            inp = convert_input_tjlf(T, input_tjlf)
            inp.RLNS[2] = x
            QL = TJLF.run_tjlf(inp)
            return TJLF.Qe(QL)
        end

        x0 = input_tjlf.RLNS[2]

        h = 1e-5
        dQe_fd = (tjlf_Qe_rlns(x0 + h) - tjlf_Qe_rlns(x0 - h)) / (2h)
        dQe_ad = ForwardDiff.derivative(tjlf_Qe_rlns, x0)

        @test isfinite(dQe_ad)
        @test isfinite(dQe_fd)
        relerr = abs(dQe_fd) > 0 ? abs(dQe_ad - dQe_fd) / abs(dQe_fd) : abs(dQe_ad - dQe_fd)
        @test relerr < 1e-4
    end

    # --- Test 4: gradient (multiple inputs) ---
    @testset "∇Qe w.r.t. [RLTS_1, RLTS_2]" begin
        function tjlf_Qe_vec(x::AbstractVector{T}) where {T<:Real}
            inp = convert_input_tjlf(T, input_tjlf)
            inp.RLTS[1] = x[1]
            inp.RLTS[2] = x[2]
            QL = TJLF.run_tjlf(inp)
            return TJLF.Qe(QL)
        end

        x0 = [input_tjlf.RLTS[1], input_tjlf.RLTS[2]]

        grad_ad = ForwardDiff.gradient(tjlf_Qe_vec, x0)

        # FD reference for each component
        h = 1e-5
        grad_fd = similar(x0)
        for i in eachindex(x0)
            xp = copy(x0); xp[i] += h
            xm = copy(x0); xm[i] -= h
            grad_fd[i] = (tjlf_Qe_vec(xp) - tjlf_Qe_vec(xm)) / (2h)
        end

        for i in eachindex(x0)
            @test isfinite(grad_ad[i])
            relerr = abs(grad_fd[i]) > 0 ? abs(grad_ad[i] - grad_fd[i]) / abs(grad_fd[i]) : abs(grad_ad[i] - grad_fd[i])
            @test relerr < 1e-4
        end
    end

    # --- Tests 5-7: SAT_RULE=1,2,3 (reuse tglf01 input, override SAT_RULE) ---
    # The IFT eigen dispatch is already compiled from SAT_RULE=0 above;
    # SAT_RULE only affects the saturation formula, so these compile quickly.
    for sat in [1, 2, 3]
        @testset "∂Qe/∂RLTS_2 (SAT_RULE=$sat)" begin
            input_tjlf_sat = deepcopy(input_tjlf)
            input_tjlf_sat.SAT_RULE = sat

            function tjlf_Qe_sat(x::T) where {T<:Real}
                inp = convert_input_tjlf(T, input_tjlf_sat)
                inp.RLTS[2] = x
                QL = TJLF.run_tjlf(inp)
                return TJLF.Qe(QL)
            end

            x0 = input_tjlf_sat.RLTS[2]
            h = 1e-5
            dQe_fd = (tjlf_Qe_sat(x0 + h) - tjlf_Qe_sat(x0 - h)) / (2h)
            dQe_ad = ForwardDiff.derivative(tjlf_Qe_sat, x0)

            @test isfinite(dQe_ad)
            @test isfinite(dQe_fd)
            relerr = abs(dQe_fd) > 0 ? abs(dQe_ad - dQe_fd) / abs(dQe_fd) : abs(dQe_ad - dQe_fd)
            @test relerr < 1e-4
        end
    end
end

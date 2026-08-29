using Test
using TJLF

# readInput validation: unknown keys and missing/invalid switches are hard errors,
# never a silent fallback (Fortran TGLF ignores unrecognized input.tglf content and
# runs its default SAT0 template — TJLF deliberately refuses instead).
@testset "readInput validation" begin
    basefile = joinpath(@__DIR__, "tglf_regression", "tglf01", "input.tglf")
    basetext = read(basefile, String)
    tmp = mktempdir()

    function read_mutated(mutation::String; drop::Union{Nothing,Regex}=nothing)
        txt = drop === nothing ? basetext : join(filter(l -> !occursin(drop, l), split(basetext, '\n')), "\n")
        f = joinpath(tmp, "input.tglf")
        write(f, txt * "\n" * mutation * "\n")
        return TJLF.readInput(f)
    end

    @testset "legacy files still parse" begin
        root = joinpath(@__DIR__, "tglf_regression")
        for d in sort(readdir(root))
            f = joinpath(root, d, "input.tglf")
            isfile(f) || continue
            @test TJLF.readInput(f) isa TJLF.InputTJLF{Float64}
        end
    end

    @testset "unknown keys are rejected" begin
        @test_throws ArgumentError read_mutated("SAT_RUL=2")                # scalar typo
        @test_throws ArgumentError read_mutated("RLST_2=3.0")               # species typo
        @test_throws ArgumentError read_mutated("MY_SPECTRUM=[1.0, 2.0]")   # vector typo
    end

    @testset "missing or invalid switches are rejected" begin
        @test_throws AssertionError read_mutated(""; drop=r"^\s*SAT_RULE\s*=")  # no SAT_RULE
        @test_throws AssertionError read_mutated("SAT_RULE=7")                  # bad domain
        @test_throws ArgumentError read_mutated("NKY=twelve")                   # unparseable value
        @test_throws AssertionError read_mutated("UNITS='GYR0'")                # UNITS typo
    end

    @testset "UNITS defaults to GYRO (Fortran default); SAT2/3 force CGYRO" begin
        @test read_mutated(""; drop=r"^\s*UNITS\s*=").UNITS == "GYRO"                # SAT0, no UNITS line
        @test read_mutated("SAT_RULE=2"; drop=r"^\s*UNITS\s*=").UNITS == "CGYRO"     # SAT2 auto-CGYRO
        # programmatic SAT2 with (default) GYRO units must be rejected, not run wrong
        inp = TJLF.readInput(basefile)
        inp.SAT_RULE = 2
        @test_throws AssertionError TJLF.checkInput(inp)
    end

    @testset "save/readInput round trip" begin
        inp = TJLF.readInput(basefile)
        f = joinpath(tmp, "roundtrip.tglf")
        TJLF.save(inp, f)
        inp2 = TJLF.readInput(f)
        for fn in fieldnames(TJLF.InputTJLF)
            startswith(String(fn), "_") && continue
            fn in (:KY_SPECTRUM, :EIGEN_SPECTRUM) && continue  # NaN-initialized, skipped on read
            @test isequal(getfield(inp2, fn), getfield(inp, fn))
        end
    end
end

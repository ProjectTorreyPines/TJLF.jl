# Plan item 1a (docs/plans/ad_typing_sysimage_container.md): measure the cost of a
# Dual solve relative to a Float64 solve, at several ForwardDiff chunk sizes N, through
# both the scalar and the Vector entry points of `run_tjlf`.
#
# Target: ratio ≈ 2 at small N (the IFT eigensolve runs once in Float64 LAPACK; the
# extra factorizations set a concrete ~2x floor), growing slowly and linearly with N
# (the O(N) derivative-propagation loops). A ratio well above that line, or one that
# grows super-linearly, indicates boxing (Union fields / dynamic-Symbol copies) or a
# fall-through to generic `eigen` on Dual inputs.
#
# Usage (from the TJLF repo root, in an env with TJLF dev'ed + ForwardDiff):
#   julia --project=<benchenv> -t 8 benchmark/ad_ratio.jl [input.tglf]
#
# Reports, for Float64 and for each N: first-call time (JIT), best-of-reps steady
# time, and steady allocations, for run_tjlf(inp) and run_tjlf([inp1..inp4]).

using TJLF
using ForwardDiff
using Printf

const INPUT_FILE = isempty(ARGS) ? joinpath(pkgdir(TJLF), "precompile", "sample_input.tglf") :
                   ARGS[1]
const CHUNK_SIZES = [1, 2, 7, 11]
const NVEC = 4      # length of the Vector entry-point batch (FUSE uses nradii ~ 7-16)
const REPS = 5

# Same shape as TJLFForwardDiffExt._convert_input_tjlf / runtests_ad.jl
function convert_input(::Type{T}, base::TJLF.InputTJLF{Float64}) where {T<:Real}
    inp = TJLF.InputTJLF{T}(base.NS, length(base.KY_SPECTRUM))
    for fn in fieldnames(TJLF.InputTJLF)
        v = getfield(base, fn)
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

function make_dual_input(base::TJLF.InputTJLF{Float64}, ::Val{N}) where {N}
    D = ForwardDiff.Dual{:adbench, Float64, N}
    inp = convert_input(D, base)
    # seed one partial so the derivative path is exercised with nonzero partials
    inp.RLTS[2] = ForwardDiff.Dual{:adbench}(base.RLTS[2],
                                             ntuple(k -> k == 1 ? 1.0 : 0.0, N))
    return inp
end

# returns (first_call_s, best_steady_s, steady_alloc_bytes)
function bench(f)
    tfirst = @elapsed f()
    best = Inf
    bytes = 0
    for _ in 1:REPS
        GC.gc()
        stats = @timed f()
        if stats.time < best
            best = stats.time
            bytes = stats.bytes
        end
    end
    return tfirst, best, bytes
end

fmt_b(b) = b > 1e9 ? @sprintf("%.2f GB", b / 1e9) : @sprintf("%.1f MB", b / 1e6)

function main()
    println("threads = $(Threads.nthreads()), input = $INPUT_FILE")
    base = TJLF.readInput(INPUT_FILE)

    rows = NamedTuple[]

    # Float64 reference
    inp64 = convert_input(Float64, base)
    f1, s1, a1 = bench(() -> TJLF.run_tjlf(inp64))
    vec64 = [convert_input(Float64, base) for _ in 1:NVEC]
    fv1, sv1, av1 = bench(() -> TJLF.run_tjlf(vec64))
    push!(rows, (; label="Float64", scalar=(f1, s1, a1), vector=(fv1, sv1, av1)))

    for N in CHUNK_SIZES
        inpD = make_dual_input(base, Val(N))
        fs, ss, as = bench(() -> TJLF.run_tjlf(inpD))
        vecD = [make_dual_input(base, Val(N)) for _ in 1:NVEC]
        fv, sv, av = bench(() -> TJLF.run_tjlf(vecD))
        push!(rows, (; label="Dual N=$N", scalar=(fs, ss, as), vector=(fv, sv, av)))
        println("  done N=$N")
    end

    s64 = rows[1].scalar[2]
    v64 = rows[1].vector[2]
    println()
    @printf("%-10s | %10s %10s %8s %10s | %10s %10s %8s %10s\n",
            "case", "sc 1st(s)", "sc best(s)", "sc ratio", "sc alloc",
            "vec 1st(s)", "vec best(s)", "vec ratio", "vec alloc")
    for r in rows
        @printf("%-10s | %10.2f %10.4f %8.2f %10s | %10.2f %10.4f %8.2f %10s\n",
                r.label,
                r.scalar[1], r.scalar[2], r.scalar[2] / s64, fmt_b(r.scalar[3]),
                r.vector[1], r.vector[2], r.vector[2] / v64, fmt_b(r.vector[3]))
    end
end

main()

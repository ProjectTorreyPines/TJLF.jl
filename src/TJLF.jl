module TJLF

using LinearAlgebra
import LinearAlgebra.LAPACK.gesv!
import LinearAlgebra.LAPACK.geev!
import LinearAlgebra.LAPACK.ggev!
using SparseArrays
using StaticArrays
using FastGaussQuadrature
using LinearMaps
import KrylovKit
import PrecompileTools

#  populated by TJLFCUDAExt.__init__() when CUDA is loaded
const _CUDA_FUNCTIONAL = Ref{Any}(() -> false)
const _CUDA_DEVICE_COUNT = Ref{Any}(() -> 0)
const _CUDA_SOLVE = Ref{Any}(nothing)
const _CUDA_LU_SOLVE = Ref{Any}(nothing)
const _CUDA_SOLVE_EIG_GRAD = Ref{Any}(nothing)

# wrappers so call sites stay the same
_cuda_functional() = _CUDA_FUNCTIONAL[]()
_cuda_device_count() = _CUDA_DEVICE_COUNT[]()
# Returns Vector{ComplexF64} of eigenvalues of B⁻¹A, computed entirely on GPU
# (CUSOLVER getrf/getrs for the solve, then Xgeev for the diagonalisation).
function _gpu_solve!(A::Matrix{ComplexF64}, B::Matrix{ComplexF64})
    _CUDA_SOLVE[] === nothing && error("CUDA extension not loaded")
    return _CUDA_SOLVE[](A, B)
end

# Solve Z x = b on the GPU (CUSOLVER getrf/getrs) and return x::Vector{ComplexF64}.
# Used for the eigenvector inverse-iteration step (Z = amat - σ·bmat); Z is overwritten.
function _gpu_lu_solve!(Z::Matrix{ComplexF64}, b::Vector{ComplexF64})
    _CUDA_LU_SOLVE[] === nothing && error("CUDA extension not loaded")
    return _CUDA_LU_SOLVE[](Z, b)
end

# GPU eigenvalues + implicit-function-theorem derivatives for the standard problem
# M = B⁻¹A (ForwardDiff.Dual support). Given the value matrices `A,B` and the per-partial
# matrices `dA[k]=∂A/∂xₖ`, `dB[k]=∂B/∂xₖ` (all host ComplexF64), this mirrors the CPU Dual
# kernel entirely on the GPU (CUSOLVER getrf/getrs for M and ∂M, Xgeev('N','V') for the
# eigenpairs, then ∂λᵢ = (R⁻¹ ∂M R)[i,i]). Returns
# `(W::Vector{ComplexF64}, dλ_re::Matrix{Float64}(n,np), dλ_im::Matrix{Float64}(n,np))`.
# The IFT diagonal R⁻¹∂M R is invariant to eigenvector column scaling, so CUSOLVER-vs-LAPACK
# normalization differences are harmless.
function _gpu_solve_eig_grad!(A::Matrix{ComplexF64}, B::Matrix{ComplexF64},
                              dA::Vector{Matrix{ComplexF64}}, dB::Vector{Matrix{ComplexF64}})
    _CUDA_SOLVE_EIG_GRAD[] === nothing && error("CUDA extension not loaded")
    return _CUDA_SOLVE_EIG_GRAD[](A, B, dA, dB)
end

# @show BLAS.get_config()
"""
    with_blas_threads(f, n::Integer)

Run `f()` with the BLAS thread count temporarily set to `n`, restoring the
previous value afterwards (even on error). No-op if BLAS is already at `n`.

TJLF parallelises over ky with Julia threads, each doing small dense LAPACK
solves; with multithreaded BLAS the kernels oversubscribe cores (e.g. 8 Julia
threads × 128 BLAS threads) and per-solve time balloons ~100x. This is applied
*scoped* around the threaded compute drivers rather than as a global default, so
loading TJLF does not change BLAS threading for the rest of a host application
(e.g. FUSE). The `n0 == n` short-circuit makes nested scopes safe: an inner
driver running inside an already-set region never touches the global BLAS state,
avoiding concurrent `set_num_threads` calls from worker threads.
"""
@inline function with_blas_threads(f, n::Integer)
    n0 = BLAS.get_num_threads()
    n0 == n && return f()
    BLAS.set_num_threads(n)
    try
        return f()
    finally
        BLAS.set_num_threads(n0)
    end
end

include("tjlf_modules.jl")
include("tjlf_read_input.jl")
include("tjlf_hermite.jl")
include("tjlf_multiscale_spectrum.jl")
include("tjlf_geometry.jl")
include("tjlf_kygrid.jl")
include("tjlf_get_uv.jl")
include("tjlf_eigensolver.jl")
include("tjlf_FLR_modules.jl")
include("tjlf_finiteLarmorRadius.jl")
include("tjlf_matrix.jl")
include("tjlf_LINEAR_SOLUTION.jl")
include("tjlf_max.jl")
include("tjlf_TRANSPORT_MODEL.jl")
include("run_tjlf.jl")

export readInput, run, get_wavefunction, pick_device
export gauss_hermite, get_sat_params, get_ky_spectrum, get_ky_spectrum_size, tjlf_TM
export sum_ky_spectrum, xgrid_functions_geo

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

# The quasi-linear spectral solve (`run_tjlf`) is by far the largest first-use JIT cost for
# downstream consumers (in the FUS3 stellarator suite it dominates a ~125 s test item, ~80% of
# which is compilation). Exercising the full read→solve path here, during precompilation, bakes
# those specializations into `TJLF.ji` so every consumer — and every parallel test worker —
# pays the compile once (cached) instead of on first call. We run both the scalar and the
# `Vector` entry points (the latter is what integrated drivers call). PrecompileTools catches
# any workload error, so this can never break precompilation.
PrecompileTools.@compile_workload begin
    input = readInput(joinpath(@__DIR__, "..", "precompile", "sample_input.tglf"))
    run_tjlf(input)
    run_tjlf([input])
end

end
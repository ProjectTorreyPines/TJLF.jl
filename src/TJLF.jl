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

#  populated by TJLFCUDAExt or TJLFEPCUDAExt.__init__() when CUDA is loaded
const _CUDA_FUNCTIONAL = Ref{Any}(() -> false)
const _CUDA_DEVICE_COUNT = Ref{Any}(() -> 0)
const _CUDA_SOLVE = Ref{Any}(nothing)

# wrappers so call sites stay the same
_cuda_functional() = _CUDA_FUNCTIONAL[]()
_cuda_device_count() = _CUDA_DEVICE_COUNT[]()
function _gpu_solve!(A::Matrix{ComplexF64}, B::Matrix{ComplexF64})
    _CUDA_SOLVE[] === nothing && error("CUDA extension not loaded")
    return _CUDA_SOLVE[](A, B)
end

# @show BLAS.get_config()
BLAS.set_num_threads(1)

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

end
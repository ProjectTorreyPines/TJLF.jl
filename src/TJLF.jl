module TJLF

using Base.Threads
using LinearAlgebra
import LinearAlgebra.LAPACK.gesv!
import LinearAlgebra.LAPACK.geev!
using SparseArrays
using StaticArrays
using FastGaussQuadrature
using LinearMaps
using KrylovKit

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

export readInput
export gauss_hermite, get_sat_params, get_ky_spectrum, get_ky_spectrum_size, tjlf_TM
export sum_ky_spectrum, xgrid_functions_geo

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end
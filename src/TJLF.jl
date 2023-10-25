# calls the fortran code as a shared library
# ccall((:main, "./src/Fortran/tglf.so"), Cvoid, () ,)

module TJLF

using LinearAlgebra
import LinearAlgebra.LAPACK: gesv!, syev!

include("tjlf_modules.jl")
include("tjlf_read_input.jl")
include("tjlf_hermite.jl")
include("intensity_sat_rules.jl")
include("tjlf_multiscale_spectrum.jl")
include("tjlf_geometry.jl")
include("tjlf_kygrid.jl")
include("tjlf_get_uv.jl")
include("tjlf_eigensolver.jl")
include("tjlf_FLR_modules.jl")
include("tjlf_finiteLarmorRadius.jl")
include("tjlf_matrix.jl")
include("tjlf_LINEAR_SOLUTION.jl")
include("tjlf_max2.jl")
include("tjlf_transport_model.jl")

export readInput
export gauss_hermite, get_sat_params, get_ky_spectrum, tjlf_TM

end
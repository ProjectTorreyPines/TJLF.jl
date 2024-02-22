module TJLFEP

using Base.Threads
using LinearAlgebra # Don't really need this for tjlf-ep
using SparseArrays
using StaticArrays

include("tjlfep_modules.jl")
include("tjlfep_read_inputs.jl")
include("EXPROconst.jl")
include("tjlfep_ky.jl")
include("tjlfep_kwscale_scan.jl")
include("mainsub.jl")
include("run_tjlfep.jl")


export profile, InputTJLF, InputTJLFEP
export readMTGLF, readTGLFEP, TJLF_map, readEXPRO
export convert_input

end 
module TJLFEP
using Revise
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
include("conv_input.jl")
include("tjlfep_complete_output.jl")
include("run_tjlfep.jl")

include("tjlfep_generate_input.jl")

#sgould this line just be the big struct?
export InputTJLFEP, profile, Options, InputTJLF
export readMTGLF, readTGLFEP, TJLF_map, readEXPRO
export convert_input, revert_input
export tjlfep_complete_output
export runTHD, runTHDs

export TJLFEP_generate_input, readline_values

end 
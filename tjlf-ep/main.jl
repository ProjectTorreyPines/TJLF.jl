using Revise

using Pkg
# Pkg.add("Revise")
# using Revise
# Pkg.add("Plots")
# Pkg.add("StaticArrays")
Pkg.activate("..")

include("TJLFEP.jl")
using .TJLFEP
using .TJLFEP: convert_input
using .TJLFEP: revert_input
include("../src/TJLF.jl")
using .TJLF
using Base.Threads
using LinearAlgebra
using Dates
BLAS.set_num_threads(1)
begin
    homedirectory = pwd()

    tglfepfilepath = homedirectory*"/../outputs/tjlfeptests/isEP3v6/input.TGLFEP"
    mtglffilepath = homedirectory*"/../outputs/tjlfeptests/isEP3v6/input.MTGLF"
    exprofilepath = homedirectory*"/../outputs/tjlfeptests/isEP3v6/input.EXPRO"

    # This now works after a brief debugging:
    runTHD(tglfepfilepath, mtglffilepath, exprofilepath)
end
# I will now run runTHDs on two examples:
#=
tglfepfilepath1 = homedirectory*"/outputs/tjlfeptests/isEP3v6/input.TGLFEP"
tglfepfilepath2 = homedirectory*"/outputs/tjlfeptests/isEP3v7/input.TGLFEP"

mtglffilepath1 = homedirectory*"/outputs/tjlfeptests/isEP3v6/input.MTGLF"
mtglffilepath2 = homedirectory*"/outputs/tjlfeptests/isEP3v7/input.MTGLF"

exprofilepath1 = homedirectory*"/outputs/tjlfeptests/isEP3v6/input.EXPRO"
exprofilepath2 = homedirectory*"/outputs/tjlfeptests/isEP3v7/input.EXPRO"

tglfeps = [tglfepfilepath1, tglfepfilepath2]
mtglfs = [mtglffilepath1, mtglffilepath2]
expros = [exprofilepath1, exprofilepath2]
=#
#widths, kymark_outs, SFmins, dpdr_crit_outs, dndr_crit_outs = runTHDs(tglfeps, mtglfs, expros)

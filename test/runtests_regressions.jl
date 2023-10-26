using Test
using Base.Filesystem
include("../src/TJLF.jl")
using ..TJLF


# saturation rule test
directory = "../outputs/tglf_regression/"
tests = readdir(directory)
for dir_name in tests
    if dir_name == ".DS_Store" continue end
    baseDirectory = satRuleDirectory*dir_name*"/"

    #******************************************************************************#************************
    # Read input.tglf
    #******************************************************************************#************************
    
    inputTJLF = readInput(baseDirectory)

    #*******************************************************************************************************
    #   start running stuff
    #*******************************************************************************************************

end
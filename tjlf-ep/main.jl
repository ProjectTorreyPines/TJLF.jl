using Revise
include("../src/TJLF.jl")
using .TJLF
using Plots
using Base.Threads
Threads.nthreads()
using LinearAlgebra
BLAS.set_num_threads(1)

include("../tjlf-ep/TJLFEP.jl")
using .TJLFEP

filename = "/Users/benagnew/gacode_add/TGLF-EP/examples/qmin1/fixw_1.88/input.TGLFEP"
inputTGLFEP = readTGLFEP(filename)

#baseDirectory = "/home/agnewb/.julia/dev/TJLF.jl/outputs/tglf_regression/tglf01/"
baseDirectory = "/Users/benagnew/TJLF.jl/outputs/tglf_regression/tglf01/"
inputTJLF = readInput(baseDirectory*"input.tglf")

outputHermite = gauss_hermite(inputTJLF)
satParams = get_sat_params(inputTJLF)
inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
QL_weight, eigenvalue, field_weight_out = tjlf_TM(inputTJLF, satParams, outputHermite)

QL_flux_out, flux_out = sum_ky_spectrum(inputTJLF, satParams, eigenvalue[:,:,1], QL_weight)

final_flux = sum(QL_flux_out,dims=1)[1,:,:]

nmodes = inputTJLF.NMODES;
fileDirectory = baseDirectory * "out.tglf.eigenvalue_spectrum";
lines = readlines(fileDirectory);
lines = split(join(lines[3:length(lines)]));
lines = [parse(Float64, l) for l in lines];

gamma = []
freq = []
for k in 1:nmodes
    push!(gamma, lines[2k-1:2*nmodes:end])
    push!(freq, lines[2k:2*nmodes:end])
end
gammaJulia = eigenvalue[:,:,1]
freqJulia = eigenvalue[:,:,2]

plot(inputTJLF.KY_SPECTRUM, freq[1]; label="Fortran")
display(plot!(inputTJLF.KY_SPECTRUM, freqJulia[1,:]; label="Julia", title="Frequency vs ky for F and J", linestyle=:dash))

plot(inputTJLF.KY_SPECTRUM, gamma[1]; label="Fortran")
display(plot!(inputTJLF.KY_SPECTRUM, gammaJulia[1,:]; label="Julia", title="Gamma vs ky for F and J", linestyle=:dash))

# The driver.jl file is not really needed and I will be deleting it soon. For now, I will be making the driver in this portion
# of the main.jl. It should help get everything in order for the most part. There are a few things the driver does that I need
# to perform (especially MPI processes), and making a separate function doesn't make too much sense to be honest.

#=============================================================================================================================#

begin
    using MPI
    include("TJLFEP.jl")
    using .TJLFEP
    using .TJLFEP: convert_input
    include("../src/TJLF.jl")
    using .TJLF

    # EXPRO definitions:

    begin
    ni = TJLFEP.exproConst.ni
    Ti = TJLFEP.exproConst.Ti
    dlnnidr = TJLFEP.exproConst.dlnnidr
    dlntidr = TJLFEP.exproConst.dlntidr
    end

    MPI.Init()
    TGLFEP_COMM_WORLD = MPI.COMM_WORLD # This is part of one of the modules in TGLFEP. Not sure what to do about that

    id = MPI.Comm_rank(MPI.COMM_WORLD)
    np = MPI.Comm_size(MPI.COMM_WORLD)

    l_print = false
    if (id == 0) 
        l_print = true 
    end

    inputEPfile = "/Users/benagnew/TJLF.jl/outputs/tglfep_tests/input.TGLFEP"
    inputMPfile = "/Users/benagnew/TJLF.jl/outputs/tglfep_tests/input.MTGLF"

    prof = TJLFEP.readMTGLF(inputMPfile)
    inputMTGLF = prof[1]
    ir_exp = prof[2]

    begin
    inputTGLFEP = TJLFEP.readTGLFEP(inputEPfile, ir_exp)

    dpdr_EP = fill(NaN, inputMTGLF.NR)
    if (inputTGLFEP.INPUT_PROFILE_METHOD == 2)
        # Got EXPRO! yay. Well the constants for IS_EP == 2 and 3 not the others though,
        for i in eachindex(dpdr_EP)
            dpdr_EP[i] = ni[inputTGLFEP.IS_EP][i]*Ti[inputTGLFEP.IS_EP][i]*(dlnnidr[inputTGLFEP.IS_EP][i]+dlntidr[inputTGLFEP.IS_EP][i])
        end
        dpdr_EP_abs = abs.(dpdr_EP)
        dpdr_EP_max = maximum(dpdr_EP_abs)
        dpdr_EP_max_loc = argmax(dpdr_EP_abs)
        n_at_max = ni[inputTGLFEP.IS_EP][dpdr_EP_max_loc]
        if (inputTGLFEP.PROCESS_IN == 5)
            for ir = 1:inputTGLFEP.SCAN_N
                println(inputTGLFEP.SCAN_N)
                inputTGLFEP.FACTOR = inputTGLFEP.FACTOR*dpdr_EP_max/dpdr_EP_abs[ir_exp[ir]]
                println(inputTGLFEP.FACTOR)
            end
        end
        inputTGLFEP.FACTOR_MAX_PROFILE .= inputTGLFEP.FACTOR
    end
    end

    # Run Mainsub:
    key = Int(id / (inputTGLFEP.SCAN_N)) 
    color = Int(id - key*(inputTGLFEP.SCAN_N))
    TJLFEP_COMM_IN = MPI.Comm_split(MPI.COMM_WORLD, color, key)

    if (inputTGLFEP.INPUT_PROFILE_METHOD == 2)
        inputTGLFEP.IR = inputTGLFEP.IR_EXP[color+1]
    else
        inputTGLFEP.IR = color + inputMTGLF.IRS
    end

    # next is the str_r which I'm skipping for now and will come back to when I have results I want to print

    inputTGLFEP.FACTOR_IN = inputTGLFEP.FACTOR[color+1]
end
# I will skip over the width statement as I believe right now that it is redundant and not really needed in Julia.
#using .TJLFEP: mainsub
#end
# This is where mainsub is now called. I will make mainsub for in case I want to work in other process-in's later
TJLFEP.mainsub(inputTGLFEP, inputMTGLF, TJLFEP_COMM_IN)

# Checking out some info:
@profview TJLFEP.mainsub(inputTGLFEP, inputMTGLF, TJLFEP_COMM_IN)
using BenchmarkTools
@btime TJLFEP.mainsub(inputTGLFEP, inputMTGLF, TJLFEP_COMM_IN)
# The following list is for unpopulated fields and things I want to look into today:
# inputTJLF.USE_AVE_ION_GRID
# inputTJLF.WIDTH_SPECTRUM
# inputTJLF.FIND_EIGEN
# inputTJLF.RLNP_CUTOFF
# inputTJLF.BETA_LOC
# inputTJLF.DAMP_PSI
# inputTJLF.DAMP_SIG
# inputTJLF.WDIA_TRAPPED

# Before I continue, I am going to make sure that the width_in terms are all correctly allotted. Optimization right now is not really
# what I am focused on. That will come later when I am closer to completion. Right now, one pass takes 265440 allocations, which is
# pretty high. TJLF itself takes up 255.41 MiB, and this is just one pass of TJLFEP running TJLF. Not even the rest of the allocations
# I need to perform on the data taken from TJLF. Speed-wise, this is, of course, why MPI is used. Allocations, however will need to probably
# be reduced. With the current trend over the 4 rounds over 500 with 1 process, this will take 66 minutes and have half a billion allocations...
# and I'm not even finished yet. This is obviously going to need to be reduced as much as possibly. Using previously-defined structs so it's
# not a pure allocation each time or using pre-definitions (There will likely be other methods, just these are the basic ones).
# Memory-wise, this is 172 GiB. That's way too large.


#= run_TJLF isolated:
begin
    using Revise
    include("../src/TJLF.jl")
    using .TJLF
    using Plots
    using Base.Threads
    Threads.nthreads()
    using LinearAlgebra
    BLAS.set_num_threads(1)
end

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
=#

# The driver.jl file is not really needed and I will be deleting it soon. For now, I will be making the driver in this portion
# of the main.jl. It should help get everything in order for the most part. There are a few things the driver does that I need
# to perform (especially MPI processes), and making a separate function doesn't make too much sense to be honest.

#=============================================================================================================================#

# In order to run TJLF-EP, MPI is used for parallel computing. Eventually I will
# create a batch method for running it, but until then, you can run it by running
# mpiexec -n # julia --project /Path/To/main.jl
# where the # is replaced by the number of processes desired 
# (generally should be >= SCAN_N)

begin
    using MPI
    include("TJLFEP.jl")
    using .TJLFEP
    using .TJLFEP: convert_input
    include("../src/TJLF.jl")
    using .TJLF

    # EXPRO definitions:

    #begin
    ni = TJLFEP.exproConst.ni
    Ti = TJLFEP.exproConst.Ti
    dlnnidr = TJLFEP.exproConst.dlnnidr
    dlntidr = TJLFEP.exproConst.dlntidr
    cs = TJLFEP.exproConst.cs
    rmin_ex = TJLFEP.exproConst.rmin_ex
    #end

    MPI.Init()
    TJLFEP_COMM_WORLD = MPI.COMM_WORLD

    id = MPI.Comm_rank(TJLFEP_COMM_WORLD)
    np = MPI.Comm_size(TJLFEP_COMM_WORLD)

    l_print = false
    if (id == 0) 
        l_print = true 
    end

    # These should be set from the working directory, but these test cases are good for now:
    inputEPfile = "/Users/benagnew/TJLF.jl/outputs/tglfep_tests/input.TGLFEP"
    inputMPfile = "/Users/benagnew/TJLF.jl/outputs/tglfep_tests/input.MTGLF"


    prof = TJLFEP.readMTGLF(inputMPfile)
    inputMTGLF = prof[1]
    ir_exp = prof[2]

    #begin
    inputTGLFEP = TJLFEP.readTGLFEP(inputEPfile, ir_exp)
    #inputTGLFEP.F_REAL
    dpdr_EP = fill(NaN, inputMTGLF.NR)
    if (inputTGLFEP.INPUT_PROFILE_METHOD == 2)
        @assert inputTGLFEP.IS_EP != 2 || inputTGLFEP.IS_EP != 3 "IS_EP value not supported in EXPRO. See EXPROconst.jl"
        for i in eachindex(dpdr_EP)
            dpdr_EP[i] = ni[inputTGLFEP.IS_EP][i]*Ti[inputTGLFEP.IS_EP][i]*(dlnnidr[inputTGLFEP.IS_EP][i]+dlntidr[inputTGLFEP.IS_EP][i])
        end
        #println(inputTGLFEP.FACTOR)
        dpdr_EP_abs = abs.(dpdr_EP)
        dpdr_EP_max = maximum(dpdr_EP_abs)
        dpdr_EP_max_loc = argmax(dpdr_EP_abs)
        n_at_max = ni[inputTGLFEP.IS_EP][dpdr_EP_max_loc]
        if (inputTGLFEP.PROCESS_IN != 5)
            for ir = 1:inputTGLFEP.SCAN_N
                #println(inputTGLFEP.SCAN_N)
                inputTGLFEP.FACTOR = inputTGLFEP.FACTOR*dpdr_EP_max/dpdr_EP_abs[ir_exp[ir]]
                #println(inputTGLFEP.FACTOR)
            end
        end
        inputTGLFEP.FACTOR_MAX_PROFILE .= inputTGLFEP.FACTOR
    end
    #end

    # Run Mainsub Setup:
    # negative color?: 2, 8, 5
    #id = 3
    #inputTGLFEP.SCAN_N = 3
    key = id / (inputTGLFEP.SCAN_N)
    key = Int(floor(key))
    color = id - key*(inputTGLFEP.SCAN_N)
    color = Int(floor(color))
    
    
    TJLFEP_COMM_IN = MPI.Comm_split(TJLFEP_COMM_WORLD, color, key)

    inputTGLFEP.F_REAL .= 1.0
    if (inputTGLFEP.REAL_FREQ == 1) 
        inputTGLFEP.F_REAL .= (cs[:]/(rmin_ex[inputMTGLF.NR]))/(2*pi*1.0e3)
    end

    if (inputTGLFEP.INPUT_PROFILE_METHOD == 2)
        # Allotting Ir_exp not from inputMTGLF.
        inputTGLFEP.IR_EXP = fill(NaN, inputTGLFEP.SCAN_N)
        for i = 1:inputTGLFEP.SCAN_N
            if (inputTGLFEP.SCAN_N != 1)
                jr_exp = inputMTGLF.IRS + floor((i-1)*(inputMTGLF.NR-inputMTGLF.IRS)/(inputTGLFEP.SCAN_N-1))
            else
                jr_exp = inputMTGLF.IRS
            end
            inputTGLFEP.IR_EXP[i] = jr_exp
        end
        # each MPI color group is assigned one ir_exp:
        inputTGLFEP.IR = inputTGLFEP.IR_EXP[color+1] 
        ir::Int = inputTGLFEP.IR
    else
        inputTGLFEP.IR = color + inputMTGLF.IRS
        ir::Int = inputTGLFEP.IR
    end

    str_r = string(Char(round(Int,ir/100) + UInt32('0'))) * 
        string(Char(round(Int,mod(ir, 100)/10) + UInt32('0'))) * 
        string(Char(mod(ir, 10) + UInt32('0')))

    inputTGLFEP.SUFFIX = "_r"*str_r

    inputTGLFEP.FACTOR_IN = inputTGLFEP.FACTOR[color+1]
end
#println(suffix)
#return 0
#println(inputTGLFEP.REAL_FREQ)
using Dates
if (id == 0)
    start_time = now()
end
growthrate = TJLFEP.mainsub(inputTGLFEP, inputMTGLF, TJLFEP_COMM_IN) # For PROCESS_IN == 5, there will be a guarantee for WIDTH_IN_FLAG being false,
if (id == 0)
    end_time = now()

    time_diff = end_time - start_time
    time_in_seconds = Dates.value(time_diff) / 1000
    println(time_in_seconds, " for np = ", np)
    println(growthrate)
end
#println(growthrate)
# Right now, width_in is being set directly from first entry from the TGLFEP input
# rather than some range of widths as TGLFEP does. This can be easily implemented but 
# most inputs won't need that for now.
#if (!inputTGLFEP.WIDTH_IN_FLAG)
#    if (id == 0)
#        a = 2
#    elseif (id < inputTGLFEP.SCAN_N)
#        b = 2
#    end
#end

# Checking out some info:
#@profview TJLFEP.mainsub(inputTGLFEP, inputMTGLF, TJLFEP_COMM_IN)
#using BenchmarkTools
#@btime TJLFEP.mainsub(inputTGLFEP, inputMTGLF, TJLFEP_COMM_IN)
# The following list is for unpopulated fields and things I want to look into today:
# inputTJLF.USE_AVE_ION_GRID
# inputTJLF.WIDTH_SPECTRUM
# inputTJLF.FIND_EIGEN
# inputTJLF.RLNP_CUTOFF
# inputTJLF.BETA_LOC
# inputTJLF.DAMP_PSI
# inputTJLF.DAMP_SIG
# inputTJLF.WDIA_TRAPPED
#using MPI
include("TJLFEP.jl")
using .TJLFEP
using .TJLFEP: convert_input
using .TJLFEP: revert_input
include("../src/TJLF.jl")
using .TJLF
using Base.Threads
using LinearAlgebra
BLAS.set_num_threads(1)

#===================================================================================================#
# main.jl is for testing the Threads version.   #
#===================================================================================================#

# Things needed:
# Number of ir (color). This is the main driver of the thread in driver.
# Each of these also has a few things that depend on n_ir. The ir chosen,
# the str_r, factor_in, and SUFFIX

ni = TJLFEP.exproConst.ni
Ti = TJLFEP.exproConst.Ti
dlnnidr = TJLFEP.exproConst.dlnnidr
dlntidr = TJLFEP.exproConst.dlntidr
cs = TJLFEP.exproConst.cs
rmin_ex = TJLFEP.exproConst.rmin_ex

# These should be set from the working directory, but these test cases are good for now:

homedir = pwd()

iEPexist::Bool = false
iMPexist::Bool = false

println(homedir)

iEPexist = isfile(homedir*"/input.TGLFEP")
iMPexist = isfile(homedir*"/input.MTGLF")

#@assert iEPexist != false "TGLFEP input file not found in directory"
#@assert iMPexist != false "MTGLF input file not found in directory"

inputEPfile = homedir*"/input.TGLFEP"
inputMPfile = homedir*"/input.MTGLF"

inputEPfile = "/Users/benagnew/TJLF.jl/outputs/tjlfeptests/isEP2copy/input.TGLFEP"
inputMPfile = "/Users/benagnew/TJLF.jl/outputs/tjlfeptests/isEP2copy/input.MTGLF"

prof = TJLFEP.readMTGLF(inputMPfile)
inputMTGLF = prof[1]
ir_exp = prof[2]

#begin
inputTGLFEP = TJLFEP.readTGLFEP(inputEPfile, ir_exp)
#inputTGLFEP.FACTOR_MAX_PROFILE
dpdr_EP = fill(NaN, inputMTGLF.NR)
if (inputTGLFEP.INPUT_PROFILE_METHOD == 2)
    #@assert inputTGLFEP.IS_EP != 2 || inputTGLFEP.IS_EP != 3 "IS_EP value not supported in EXPRO. See EXPROconst.jl"
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
println(dpdr_EP)

# After the split:

inputTGLFEP.F_REAL .= 1.0
if (inputTGLFEP.REAL_FREQ == 1) 
    inputTGLFEP.F_REAL .= (cs[:]/(rmin_ex[inputMTGLF.NR]))/(2*pi*1.0e3)
end

if (inputTGLFEP.INPUT_PROFILE_METHOD == 2)
    # Allotting Ir_exp not from inputMTGLF.
    inputTGLFEP.IR_EXP = fill(0, inputTGLFEP.SCAN_N)
    for i = 1:inputTGLFEP.SCAN_N
        if (inputTGLFEP.SCAN_N != 1)
            jr_exp = inputMTGLF.IRS + floor((i-1)*(inputMTGLF.NR-inputMTGLF.IRS)/(inputTGLFEP.SCAN_N-1))
        else
            jr_exp = inputMTGLF.IRS
        end
        inputTGLFEP.IR_EXP[i] = jr_exp
    end
end


# Needs to be corrected later (in mainsub):
#=
    # each MPI color group is assigned one ir_exp:
    # !correct color!:
    inputTGLFEP.IR = inputTGLFEP.IR_EXP[color+1] 
    ir::Int = inputTGLFEP.IR
    else
    inputTGLFEP.IR = color + inputMTGLF.IRS
    ir::Int = inputTGLFEP.IR

str_r = string(Char(round(Int,ir/100) + UInt32('0'))) * 
    string(Char(round(Int,mod(ir, 100)/10) + UInt32('0'))) * 
    string(Char(mod(ir, 10) + UInt32('0')))

inputTGLFEP.SUFFIX = "_r"*str_r

inputTGLFEP.FACTOR_IN = inputTGLFEP.FACTOR[color+1]
=#

# Things that need to be set for the function of the driver:
# str_r & SUFFIX
# FACTOR_IN
# color & id (sort of)

n_ir = inputTGLFEP.SCAN_N
T1 = deepcopy(inputTGLFEP)
T2 = deepcopy(T1)
T3 = deepcopy(T2)
arrTGLFEP = [T1, T2, T3]
arrMTGLF = fill(inputMTGLF, n_ir)
arrgrowth = fill(fill(NaN,(5, 10, 10, inputTGLFEP.NMODES)), n_ir)

using Dates

start_time = now()

Threads.@threads for i in 1:n_ir
    #try

        arrTGLFEP[i].IR = arrTGLFEP[i].IR_EXP[i]
        #println(arrTGLFEP[i].IR)
        ir = arrTGLFEP[i].IR
        str_r = string(Char(round(Int,ir/100) + UInt32('0'))) * 
        string(Char(round(Int,mod(ir, 100)/10) + UInt32('0'))) * 
        string(Char(mod(ir, 10) + UInt32('0')))
        arrTGLFEP[i].SUFFIX = "_r"*str_r
        #println(arrTGLFEP[i].SUFFIX)
        arrTGLFEP[i].FACTOR_IN = arrTGLFEP[i].FACTOR[i]
        input1 = arrTGLFEP[i]
        cond1 = (i > 1 && arrTGLFEP[i] == arrTGLFEP[i-1])
        cond2 = (i > 2 && arrTGLFEP[i] == arrTGLFEP[i-2])
        if (cond1 || cond2)
            if (cond1)
                println("arrTGLFEP at ", i, " = ", i-1)
            elseif (cond2)
                println("arrTGLFEP at ", i, " = ", i-2)
            end
        else
            println("All good. 1 != 2 and 2 != 3 and 3 != 1")
        end
        input2 = arrMTGLF[i]
        arrgrowth[i], arrTGLFEP[i], arrMTGLF[i] = TJLFEP.mainsub(input1, input2)
    #catch
    #end
end



end_time = now()

time_diff = end_time - start_time
time_in_seconds = Dates.value(time_diff) / 1000
println(time_in_seconds)

arrTGLFEP[3].IR
arrTGLFEP[3].FACTOR_IN
#println(growthrate)
# Right now, width_in is being set directly from first entry from the TGLFEP input
# rather than some range of widths as TGLFEP does. This can be easily implemented but 
# most inputs won't need that for now.

#println(inputTGLFEP.WIDTH_IN, " for id & color ", id_world, " ", color)

kymark_out::Vector{Float64} = fill(NaN, inputTGLFEP.SCAN_N)
width::Vector{Float64} = fill(NaN, inputTGLFEP.SCAN_N)
#local buf_width::Vector{Float64} = fill(NaN, 1)
#local buf_kymark::Vector{Float64} = fill(NaN, 1)

if (!inputTGLFEP.WIDTH_IN_FLAG)
    # Non-MPI:
    # There are only "3" processes in Threads -- 
    for i = 1:n_ir
        width[i] = arrTGLFEP[i].WIDTH_IN
        kymark_out[i] = arrTGLFEP[i].KYMARK
    end

    # Original MPI:
    #=
    #println("id & color: ", id_world, " ", color)
    MPI.Barrier(MPI.COMM_WORLD)
    println("Here! ", id_world, " ", color)
    if (id_world < inputTGLFEP.SCAN_N && id_world != 0)
        # This lines up with the Fortran.
        MPI.Send([inputTGLFEP.WIDTH_IN], 0, id_world, MPI.COMM_WORLD) # No, dest=0 is correct.
        MPI.Send([inputTGLFEP.KYMARK], 0, id_world+inputTGLFEP.SCAN_N, MPI.COMM_WORLD)
        
        println("id, ir, sent WIDTH_IN and KYMARK: ", id_world, " : ", inputTGLFEP.IR, " : ", inputTGLFEP.WIDTH_IN, " : ", inputTGLFEP.KYMARK)
        #println(buf_width, " before ", id, " & ", color)
    end
    MPI.Barrier(MPI.COMM_WORLD)
    println("Here! ", id_world, " ", color)
    if (id_world == 0)
        width[1] = inputTGLFEP.WIDTH_IN
        #println(width[1])
        kymark_out[1] = inputTGLFEP.KYMARK
        # I believe the reason Scan_n-1 works here is because all of the processes will have been reduced, and thus will not be of concern
        # This does mean that there will be things sent that are not received.
        for i = 1:inputTGLFEP.SCAN_N-1 
            buf_width = [NaN]
            #println([width[i+1]], " before ", id_world, " & ", color, " & ", i)
            buf_kymark = [NaN]
            MPI.Recv!(buf_width, i, i, MPI.COMM_WORLD)
            MPI.Recv!(buf_kymark, i, i+inputTGLFEP.SCAN_N, MPI.COMM_WORLD)
            width[i+1] = buf_width[1]
            kymark_out[i+1] = buf_kymark[1]
            #println([width[i+1]], " before ", id_world, " & ", color, " & ", i)
            #width[i+1] = buf_width[1]
            #kymark_out[i+1] = buf_kymark[1]
        end
    end
    =#
end
arrTGLFEP[1].WIDTH_IN

# As we have one process now, we don't need to check for id == 0:
inputTGLFEP = arrTGLFEP[1]
io2 = open("out.TGLFEP", "w")
println(io2, "process_in = ", inputTGLFEP.PROCESS_IN)

if (inputTGLFEP.PROCESS_IN <= 1) println(io2, "mode_in = ", inputTGLFEP.MODE_IN) end
if ((inputTGLFEP.PROCESS_IN == 4) || (inputTGLFEP.PROCESS_IN == 5)) println(io2, "threshold_flag = ", inputTGLFEP.THRESHOLD_FLAG) end

println(io2, "ky_mode = ", inputTGLFEP.KY_MODEL)
println(io2, "--------------------------------------------------------------")
println(io2, "scan_n = ", inputTGLFEP.SCAN_N)
println(io2, "irs = ", inputTGLFEP.IRS)
println(io2, "n_basis = ", inputTGLFEP.N_BASIS)
println(io2, "scan_method = ", inputTGLFEP.SCAN_METHOD)

if (inputTGLFEP.WIDTH_IN_FLAG)
    println(io2, "ir,  width")
    for i = 1:inputTGLFEP.SCAN_N
        println(io2, inputTGLFEP.IRS+i-1, " ", width[i])
    end
else
    println(io2, "ir,  width,  kymark")
    for i = 1:inputTGLFEP.SCAN_N
        println(io2, inputTGLFEP.IRS+i-1, " ", width[i], " ", kymark_out[i])
    end
end

println(io2, "--------------------------------------------------------------")
println(io2, "factor_in_profile = ", inputTGLFEP.FACTOR_IN_PROFILE)
if (inputTGLFEP.FACTOR_IN_PROFILE)
    for i = 1:inputTGLFEP.SCAN_N
        println(io2, inputTGLFEP.FACTOR[i])
    end
else
    println(io2, inputTGLFEP.FACTOR[1])
end

println(io2, "width_in_flag = ", inputTGLFEP.WIDTH_IN_FLAG)
if (!inputTGLFEP.WIDTH_IN_FLAG) println(io2, "width_min = ", inputTGLFEP.WIDTH_MIN, " width_max = ", inputTGLFEP.WIDTH_MAX) end
close(io2)


# Lastly is the printing of the density threshold:
#local buf_factor::Vector{Float64} = fill(NaN, 1)
if ((inputTGLFEP.PROCESS_IN == 4) || (inputTGLFEP.PROCESS_IN == 5))
    # Threads Process won't happen until later:


    # Original MPI:
    #=if (id_world < inputTGLFEP.SCAN_N && id_world != 0)
        if (inputTGLFEP.THRESHOLD_FLAG == 0)
            #buf_factor = [inputTGLFEP.FACTOR_IN]
            MPI.Send([inputTGLFEP.FACTOR_IN], 0, id_world, MPI.COMM_WORLD)
            #inputTGLFEP.FACTOR_IN = buf_factor[1]
        end
        #MPI.Send(inputTGLFEP) only process-In 4 here.
    end # id == 0
    MPI.Barrier(TJLFEP_COMM_IN)=#
    
    io3 = open("out.TGLFEP", "a")
    println(io3, "**************************************************************")
    println(io3, "************** The critical EP density gradient **************")
    println(io3, "**************************************************************")

    SFmin = fill(NaN, inputTGLFEP.SCAN_N)
    SFmin_out = fill(NaN, inputMTGLF.NR)
    dndr_crit = fill(NaN, inputTGLFEP.SCAN_N)
    dndr_crit_out = fill(NaN, inputMTGLF.NR)
    dpdr_crit = fill(NaN, inputTGLFEP.SCAN_N)
    dpdr_crit_out = fill(NaN, inputMTGLF.NR)

    if (inputTGLFEP.THRESHOLD_FLAG == 0)
        # Threads Correct:
        for i = 1:n_ir
            SFmin[i] = arrTGLFEP[i].FACTOR_IN
        end
        # MPI Original:
        #SFmin[1] = inputTGLFEP.FACTOR_IN
        #println(io3, inputTGLFEP.FACTOR_IN, " factor_in before")
        #println("Before MPI.Recv! for factor_in.")
        #=for i = 1:inputTGLFEP.SCAN_N-1
            buf_factor = [NaN]
            MPI.Recv!(buf_factor, i, i, MPI.COMM_WORLD)
            SFmin[i+1] = buf_factor[1]
            #println(io3, SFmin[i_1], " factor_in after and ", inputTGLFEP.FACTOR_IN, " buf_factor after")
            #println(io3, SFmin) # before buf_factor
            #SFmin[i+1] = buf_factor[1]
            #println(io3, SFmin) # after buf_factor, before comp.out
            
        end=#
        println("After MPI.Recv! for factor_in")
        println(io3, "--------------------------------------------------------------")
        println(io3, "SFmin")
        
        # Next is TGLFEP_complete_output(SFmin, SFmin_out, ir_min, ir_max, l_accept_profile)
        SFmin, SFmin_out, ir_min, ir_max, l_accept_profile = tjlfep_complete_output(SFmin, inputTGLFEP, inputMTGLF)
        println(io3, SFmin, " SFmin after buf and coutput") # after comp.out
        #println(io3, inputTGLFEP.FACTOR_MAX_PROFILE)
        if ((ir_min-inputTGLFEP.IRS+1) > 1)
            SFmin[1:ir_min-inputTGLFEP.IRS] .= inputTGLFEP.FACTOR_MAX_PROFILE[1:ir_min-inputTGLFEP.IRS]
            if (inputTGLFEP.IRS > 1) SFmin_out[1:inputTGLFEP.IRS-1] .= inputTGLFEP.FACTOR_MAX_PROFILE[1] end
            SFmin_out[inputTGLFEP.IRS:ir_min-1] .= inputTGLFEP.FACTOR_MAX_PROFILE[1:ir_min-inputTGLFEP.IRS]
        end
        if ((ir_max-inputTGLFEP.IRS+1) < inputTGLFEP.SCAN_N)
            SFmin[ir_max-inputTGLFEP.IRS+2:inputTGLFEP.SCAN_N] .= inputTGLFEP.FACTOR_MAX_PROFILE[ir_max-inputTGLFEP.IRS+2:inputTGLFEP.SCAN_N]
            if (inputTGLFEP.IRS+inputTGLFEP.SCAN_N-1 < inputMTGLF.NR) SFmin_out[inputTGLFEP.IRS+inputTGLFEP.SCAN_N:inputMTGLF.NR] .= inputTGLFEP.FACTOR_MAX_PROFILE[inputTGLFEP.SCAN_N] end
            SFmin_out[ir_max+1:inputTGLFEP.IRS+inputTGLFEP.SCAN_N-1] .= inputTGLFEP.FACTOR_MAX_PROFILE[ir_max-inputTGLFEP.IRS+2:inputTGLFEP.SCAN_N]
        end
        println(io3, SFmin, " SFmin after Max assign")
        for i = 1:inputTGLFEP.SCAN_N
            if (l_accept_profile[i])
                println(io3, SFmin[i])
            else
                println(io3, SFmin_out[i+inputTGLFEP.IRS-1], "   (*)")
            end
        end

        if (inputTGLFEP.INPUT_PROFILE_METHOD == 2)
            dndr_crit .= 10000.0
            for i = 1:inputTGLFEP.SCAN_N
                if (SFmin[i] < 9000.0)
                    dndr_crit[i] = SFmin[i]*ni[inputTGLFEP.IS_EP][Int(inputTGLFEP.IR_EXP[i])]*dlnnidr[inputTGLFEP.IS_EP][Int(inputTGLFEP.IR_EXP[i])]
                elseif ((i < ir_min-inputTGLFEP.IRS+1) || (i > ir_max-inputTGLFEP.IRS+1))
                    dndr_crit[i] = inputTGLFEP.FACTOR_MAX_PROFILE[i]*ni[inputTGLFEP.IS_EP][Int(inputTGLFEP.IR_EXP[i])]*dlnnidr[inputTGLFEP.IS_EP][Int(inputTGLFEP.IR_EXP[i])]
                end
            end
            dndr_crit, dndr_crit_out, ir_dum_1, ir_dum_2, l_accept_profile = tjlfep_complete_output(dndr_crit, inputTGLFEP, inputMTGLF)
            io4 = open("alpha_dndr_crit.input", "w")
            println(io4, "Density critical gradient (10^19/m^4)")
            println(io4, dndr_crit_out)
            close(io4)
        end

        if (inputTGLFEP.INPUT_PROFILE_METHOD == 2)
            dpdr_crit .= 10000.0
            dpdr_EP[:] .= ni[inputTGLFEP.IS_EP][:].*Ti[inputTGLFEP.IS_EP][:].*(dlnnidr[inputTGLFEP.IS_EP][:].+dlntidr[inputTGLFEP.IS_EP][:]).*0.16022
            for i = 1:inputTGLFEP.SCAN_N
                if (SFmin[i] < 9000.0)
                    if ((inputTGLFEP.PROCESS_IN == 4) || (inputTGLFEP.PROCESS_IN == 5))
                        case = inputTGLFEP.SCAN_METHOD
                        if (case == 1)
                            dpdr_scale = SFmin[i]
                        elseif (case == 2)
                            dpdr_scale = ((SFmin[i]*dlnnidr[inputTGLFEP.IS_EP][inputTGLFEP.IR_EXP[i]]+dlntidr[inputTGLFEP.IS_EP][inputTGLFEP.IR_EXP[i]]) /
                            (dlnnidr[inputTGLFEP.IS_EP][inputTGLFEP.IR_EXP[i]]+dlntidr[inputTGLFEP.IS_EP][inputTGLFEP.IR_EXP[i]]))
                        end
                        dpdr_crit[i] = dpdr_scale*dpdr_EP[inputTGLFEP.IR_EXP[i]]
                    end # 4 || 5
                end # < 9000
            end # over scan_n
            dpdr_crit, dpdr_crit_out, ir_dum_1, ir_dum_2, l_accept_profile = tjlfep_complete_output(dpdr_crit, inputTGLFEP, inputMTGLF)
            io5 = open("alpha_dpdr_crit.input", "w")
            println(io5, "Pressure critical gradient (10 kPa/m)")
            println(io5, dpdr_crit_out)
            close(io5)
        end # end prof. method 2

        println(io3, "--------------------------------------------------------------")
        println(io3, "The EP density threshold n_EP/n_e (%) for gamma_AE = 0")
        for i = 1:inputTGLFEP.SCAN_N
            println(io3, SFmin[i]*inputMTGLF.AS[inputTGLFEP.IRS+i-1, inputMTGLF.IS]*100.0) #percent
        end
        
        println(io3, "--------------------------------------------------------------")
        println(io3, "The EP beta crit (%) = beta_e*(n_EP_th/n_e)*(T_EP/T_e)")
        for i = 1:inputTGLFEP.SCAN_N
            if (inputMTGLF.GEOMETRY_FLAG == 0)
                println(io3, SFmin[i]*inputMTGLF.BETAE[inputTGLFEP.IRS+i-1]*100.0*inputMTGLF.AS[inputTGLFEP.IRS+i-1, inputMTGLF.IS]*inputMTGLF.TAUS[inputTGLFEP.IRS+i-1, inputMTGLF.IS]) #percent
            else
                println(io3, SFmin[i]*inputMTGLF.BETAE[inputTGLFEP.IRS+i-1]*100.0*inputMTGLF.AS[inputTGLFEP.IRS+i-1, inputMTGLF.IS]*inputMTGLF.TAUS[inputTGLFEP.IRS+i-1, inputMTGLF.IS]*inputMTGLF.KAPPA[inputTGLFEP.IRS+i-1]^2) #percent
            end
        end

        # there is a process_in == 4 addition I won't be doing quite yet.
    else # ThreshFlag != 0
        # Skipping for now as I want to test just threshold flag == 0 first
    end # ThreshFlag
end # process 4 || 5
close(io3)
"""
readMTGLF: Directly converts a known MTGLF file into a profile struct. It used to do this directly to an inputTJLF struct, but I am
going to follow more directly what TGLFEP does first and then come back and make alterations later if I so decide.

Inputs: filepath of input.MTGLF

Outputs: InputTJLF struct to be run through TJLF (as is done in TGLFEP)
"""
function readMTGLF(filename::String)
    open(filename)
    lines = readlines(filename)

    # This will be all very similar to TJLF for the most part. I'll check for NS and NKY (? - That's not used in TGLFEP so maybe the defualt)
    dflt = true
    ns = -1
    nr = -1
    
    for line in lines[1:100]
        line = split(line, "\n")
        line = split(line[1],"=")
        for i in 1:2
            line[i] = strip(line[i])
        end
        
        if line[1] == "NS"
            ns = parse(Int, strip(line[2]))
        elseif line[1] == "NR"
            nr = parse(Int, strip(line[2]))
        end
    end
    # make sure ns is defined
    @assert ns!=-1 "did not find NS in $filename"
    @assert nr!=-1 "din not find NR in $filename"
    
    #Translating MTGLF to usable InputTJLF form:
    inputMTGLF = profile{Float64}(nr, ns)
    irexp::Vector{Int64} = []
    
    for line in lines[1:length(lines)]
        line = split(line, "\n")
        line = split(line[1],"=")
        # Reformat for ease of access:
        for i in 1:2
            line[i] = strip(line[i])
        end

        check = (match(r"\d\d_\d", line[1]) !== nothing) || (match(r"\d\d\d_\d", line[1]) !== nothing) || (match(r"\d_\d", line[1]) !== nothing)
        vcheck = (match(r" \d", line[1]) !== nothing) || (match(r"_\d", line[1]) !== nothing)

        twoName = ["OMEGA_TAE", "RHO_STAR", "S_ZETA", "S_DELTA", "S_KAPPA", "P_PRIME", "Q_PRIME", "B_UNIT", "IR_EXP"]
        delFields = ["IR"]
        if (check == true) # Matrix values
            val = parse(Float64, line[2])
            line = split(line[1], "_")
            for i in 1:3
                line[i] = strip(line[i])
            end
            
            if (!contains(line[2], "SHEAR"))
                speciesField = Symbol(line[1])
                Index1 = line[2]
                Index2 = line[3]
            else # Any non-singular named will be here basically.
                speciesField = Symbol("VPAR_SHEAR")
                Index1 = line[3]
                Index2 = line[4]
            end
            getfield(inputMTGLF, speciesField)[parse(Int, Index1), parse(Int, Index2)] = val
            # Bam, done with that part. Next to vectors.

        elseif (check == false && vcheck == true) # Vector values
            val = parse(Float64, line[2])
            speciesField = strip(replace(replace(replace(line[1], r" \d"=>""), r"\d"=>""), r"\d"=>""), ['_', ' '])
            
            line = split(line[1], "_")
            
            if (speciesField ∈ twoName)
                speciesIndex = line[3]
            else
                speciesIndex = line[2]
            end
            if (speciesField != "IR_EXP")
                speciesField = Symbol(speciesField)
                getfield(inputMTGLF, speciesField)[parse(Int,speciesIndex)] = val
            else
                append!(irexp, Int(val))
            end
            # The ALPHA vector is only used for the input.profile method apparently. I cannot find it anywhere in TGLFEP besides
            # the method where it is listed in input.profile explicitly.
        else # Non-Vectors/Matrices

            field = line[1]
        
            if (contains(line[2], 'T') || contains(line[2], 'F'))
                val = lowercase(strip(line[2], ['\'','.'])) == "t"
            elseif (!contains(line[2], '.'))
                val = parse(Int64, line[2])
            else
                val = parse(Float64, line[2])
            end
            if (field ∈ delFields) continue
            else
                field = Symbol(field)
                setfield!(inputMTGLF, field, val)
            end
        end
    end
    # That's all for now folks!
    return inputMTGLF, irexp
end
"""
readprofile is meant for the input.profile method of insertion
"""
function readprofile(filename::String)

end

"""
readTGLFEP extracts the values needed from the input.TGLFEP file

Inputs: filename

Outputs: Options struct
"""
function readTGLFEP(filename::String, ir_exp::Vector{Int64})
    open(filename)
    lines = readlines(filename)

    # Required values to extract BEFORE assignment:
    nscan_in = -1
    widthin = -1
    ky_model = -1 # For assigning n_toroidal
    process_in = -1 # For assigning nn
    threshold_flag = -1 # For assigning nn
    
    for line in lines[1:length(lines)]

        # This would need to be adjusted for WIDTH_IN_FLAG = false, which is
        # not a covered case in TGLF-EP (but is in OMFIT?)

        if contains(line, "") && !contains(line, " ") continue end
        line = split(line, "\n")

        if contains(line[1], "SCAN_N")
            line = split(line[1])
            nscan_in = parse(Int, strip(line[1]))
        elseif contains(line[1], "WIDTH_IN_FLAG")
            line = split(line[1])
            widthin = parse(Bool, strip(line[1], '.'))
        elseif contains(line[1], "KY_MODEL")
            line = split(line[1])
            ky_model = parse(Int, strip(line[1]))
        elseif contains(line[1], "PROCESS_IN")
            line = split(line[1])
            process_in = parse(Int, strip(line[1]))
        elseif contains(line[1], "THRESHOLD_FLAG")
            line = split(line[1])
            threshold_flag = parse(Int, strip(line[1]))
        end

    end
    
    @assert nscan_in != -1 "SCAN_N not found in TGLFEP input"
    @assert widthin != -1 "WIDTH_IN_FLAG not found in TGLFEP input"
    @assert ky_model != -1 "KY_MODEL not found in TGLFEP input"
    @assert process_in != -1 "PROCESS_IN not found in TGLFEP input"
    @assert threshold_flag != -1 "THRESHOLD_FLAG not found in TGLFEP input"

    # Now, nscan_in, widthin are assigned and ready.

#=  # Assign nn:
if ((process_in == 4) || (process_in == 5))
    if (threshold_flag == 0)
        nn = 5
    else
        nn = 15
    end
end
=#
    # See TGLFEP_interface.f90:
    jtscale_max = 1
    nmodes = 4

    # NR is derived from the profile! This means that in order to read TJLFEP in Julia, we have to call the previous inputMTGLF
    # function call's NR. It will now be a required input for this function. 
    nr = 201
    nn = 5
    inputTJLFEP = Options{Float64}(nscan_in, widthin, nn, nr, jtscale_max, nmodes)
    # println("inputTJLFEP.IR_EXP = ",inputTJLFEP.IR_EXP)
    # println("ir_exp", ir_exp)
    inputTJLFEP.IR_EXP = ir_exp
    inputTJLFEP.NMODES = nmodes
   


    for line in lines[1:length(lines)]
        line = replace(line, "   "=>" ") # Don't ask lol
        line = replace(line, "  "=>" ")

        line = split(line, "\n")
        if contains(line[1], "") && !contains(line[1], " ") continue end # && !contains(line[1], ".") may be needed for if WIDTH_IN can change in the vector. I am yet to see a case of this.

        line = split(line[1], " ")
        # Now I will have all input parameters I want. A space must exist between the fields. It is better if it is just one space but can be UP TO a tab (3 spaces in VSCode where I am editing this).
        
        vecFields = ["WIDTH", "FACTOR"] # For now...
        if line[2] ∈ vecFields 
            #if ()
            field = Symbol(line[2])
            getfield(inputTJLFEP, field) .= [parse(Float64,line[1])]
        else
            field = Symbol(line[2])
            if line[1][1] == '.'
                val = lowercase(strip(line[1], ['\'','.'])) == "true"
            elseif !contains(line[1], '.')
                val = parse(Int, line[1])
            else
                val = parse(Float64, line[1])
            end

            try
                println("field: $field $(contains(string(field),"THETA_2_THRESH"))")
                if contains(string(field),"THETA_2_THRESH")
                    setfield!(inputTJLFEP, Symbol("THETA_SQ_THRESH"), val)
                else
                    setfield!(inputTJLFEP, field, val)
                end
            catch
                throw(error(field))
            end
        end
    end
    # This function is for the most part done. There should be consideration for default conditions however...
    # Something similar to the checkInput function created for TJLF but with the option to define a default or throw an error back if it's not populated. In fact, some things will have to be
    # unpopulated depending on whether WIDTH_IN_FLAG is true or false, e.g..

    # There are a few fields that need to be manually set. I will list them here:

    if (ky_model == 0)
        inputTJLFEP.NTOROIDAL = 4
    else
        inputTJLFEP.NTOROIDAL = 3
    end
        
    if (process_in == 4 || process_in == 5)
        inputTJLFEP.NN = nn
    end

    if (!inputTJLFEP.FACTOR_IN_PROFILE)
        inputTJLFEP.FACTOR = fill(inputTJLFEP.FACTOR_IN, nscan_in)
    end
    inputTJLFEP.FACTOR_MAX_PROFILE = inputTJLFEP.FACTOR

    return inputTJLFEP
end

"""
TJLF_map: This directly maps any needed inputs from the InputTJLFEP and profile structs that were obtained
from the input.MTGLF (or input.profile) and input.TGLFEP files in tjlfep_read_inputs.jl (above). This is done so
TJLF can be run very easily. The InputTJLF struct is defined above as needed in order to perform this.

inputs: InputTJLFEP from TGLFEP input file; profile from MTGLF profile file

Outputs: InputTJLF struct ready for usage in running TJLF. 
"""
#Temporary using statements:
#include("../tjlf-ep/TJLFEP.jl")
#using .TJLFEP

function TJLF_map(inputsEP::Options{Float64}, inputsPR::profile{Float64})
    # Access the fields like this:
    # inputsOptions = inputsEP.Options
    # profile = inputsEP.profile   
     #=
    # Temp Defs:
    color = 0
    kyhat_in = 3
    # Temp Struct Inputs:
    filename = "/Users/benagnew/TJLF.jl/outputs/tglfep_tests/input.MTGLF"
    temp = readMTGLF(filename)
    inputsPR = temp[1]
    irexp2 = temp[2]
    filename = "/Users/benagnew/gacode_add/sample-rundir_2/input.TGLFEP"
    inputsEP = readTGLFEP(filename, irexp2)
    inputsEP.IR = inputsEP.IR_EXP[color+1]
    inputsEP.MODE_IN = 2
    inputsEP.KY_MODEL = 3
    =#
    inputsEP.MODE_IN = 2
    inputsEP.KY_MODEL = 3

    # Okay finally I can do this lol:
    inputTJLF = InputTJLF{Float64}(inputsPR.NS, 12, true) # It is being set to the default...
    if (inputsEP.IR < 1 || inputsEP.IR > inputsPR.NR)
        println("ir isn't within range")
        return 1
    end
    inputTJLF.SIGN_BT = inputsPR.SIGN_BT
    inputTJLF.SIGN_IT = inputsPR.SIGN_IT

    inputTJLF.SAT_RULE = 0

    inputTJLF.NS = inputsPR.NS
    ns = inputsPR.NS

    is = inputsEP.IS_EP + 1
    inputsPR.IS = is
    ir = inputsEP.IR

    #TJLF deletes GEOMETRY_FLAG so this is redundant:
    #inputTJLF.GEOMETRY_FLAG = inputsEP.GEOMETRY_FLAG

    inputTJLF.ZS = inputsPR.ZS
    inputTJLF.MASS = inputsPR.MASS
    inputTJLF.AS = inputsPR.AS[ir, :] # Check read_inputs for these
    inputTJLF.TAUS = inputsPR.TAUS[ir, :]
    
    inputTJLF.ZS[1] = -1.0
    
    inputsEP.FACTOR_MAX = 0.5*1.0/(inputTJLF.ZS[is]*inputsPR.AS[ir, is])
    if (inputsEP.SCAN_METHOD == 2)
        inputsEP.FACTOR_MAX = 1.0E3
    end
    inputsEP.FACTOR_IN
    if (inputsEP.FACTOR_IN < 0)
        inputsEP.FACTOR_IN = 0
    end
    if (inputsEP.FACTOR_IN > inputsEP.FACTOR_MAX)
        inputsEP.FACTOR_IN = inputsEP.FACTOR_MAX
    end

    # Trying to correct for rounding point errors:
    # e.g. 0.1000000001 -> 0.1

    inputsEP.FACTOR_IN = round(100000*inputsEP.FACTOR_IN)/100000

    # Factor_in is used to scale only the ion after the energetic species:
    if (inputsEP.SCAN_METHOD == 1)
        inputTJLF.AS[is] = inputsPR.AS[ir, is]*inputsEP.FACTOR_IN
    end
    # I believe this is a correction for quasineutrality? AS is the ratio of density of the species to the electron

    sum0 = 0.0
    for i = 2:ns
        if (i != is) 
            sum0 = sum0 + inputTJLF.ZS[i]*inputTJLF.AS[i]
        end
    end
    if (inputsEP.IR == 2 && false)
        # println("======")
        # println(inputsPR.A_QN)
        # println(inputTJLF.ZS)
    end
    inputsPR.A_QN = (1.0 - inputTJLF.ZS[is]*inputTJLF.AS[is]) / sum0

    if (inputsEP.IR == 2 && false)
        # println(sum0)
        # println(inputsPR.A_QN)
        # println(inputTJLF.AS)
        # println("======")
    end
    for i = 2:ns
        if (i != is)
            inputTJLF.AS[i] = inputsPR.A_QN*inputTJLF.AS[i]
        end
    end
    # if (inputsEP.IR == 2 && false)
    #     println(inputTJLF.AS)
    
    #println("is, FACTOR_MAX, FACTOR_IN, sum0, AS:")
    #println(is, " ", inputsEP.FACTOR_MAX, " ", inputsEP.FACTOR_IN, " ", sum0, " ", inputTJLF.AS)

    if (inputsEP.MODE_IN == 2) # EP drive only
        for i = 1:ns
            if (i != is)
                inputTJLF.RLNS[i] = 1.0e-6
                inputTJLF.RLTS[i] = 1.0e-6
            else
                inputTJLF.RLNS[i] = inputsPR.RLNS[ir, i]
                inputTJLF.RLTS[i] = inputsPR.RLTS[ir, i]
            end
        end
    else
        for i = 1:ns
            inputTJLF.RLNS[i] = inputsPR.RLNS[ir, i]
            inputTJLF.RLTS[i] = inputsPR.RLTS[ir, i]
        end
    end

    # EP species
    # etc...

    if (inputsEP.MODE_IN == 3)
        for i = is:ns
            inputTJLF.RLNS[i] = 1.0e-5
            inputTJLF.RLTS[i] = 1.0e-5
        end
    else
        inputTJLF.RLNS[is] = inputsPR.RLNS[ir, is]*inputsEP.SCAN_FACTOR
        inputTJLF.FILTER = 0.0
    end

    if (inputsEP.SCAN_METHOD == 2) # SCAN_METHOD 1 or 2 is chosen for whether you are applying the scaling to the density gradient or the density.
        inputTJLF.RLNS[is] = inputsEP.FACTOR_IN*inputsPR.RLNS[ir, is]
    end

    # Geometry: TJLF can only recognize GEOMETRY_FLAG == 1 so I will skip geoflag == 0 for now

    if (inputsPR.GEOMETRY_FLAG == 0)
        
    end
    if (inputsPR.GEOMETRY_FLAG == 1)
        inputTJLF.RMIN_LOC = inputsPR.RMIN[ir]
        inputTJLF.RMAJ_LOC = inputsPR.RMAJ[ir]
        inputTJLF.ZMAJ_LOC = 0.0
        inputTJLF.DRMAJDX_LOC = inputsPR.SHIFT[ir]
        inputTJLF.DZMAJDX_LOC = 0.0
        inputTJLF.KAPPA_LOC = inputsPR.KAPPA[ir]
        inputTJLF.S_KAPPA_LOC = inputsPR.S_KAPPA[ir]
        inputTJLF.DELTA_LOC = inputsPR.DELTA[ir]
        inputTJLF.S_DELTA_LOC = inputsPR.S_DELTA[ir]
        inputTJLF.ZETA_LOC = inputsPR.ZETA[ir]
        inputTJLF.S_ZETA_LOC = inputsPR.S_ZETA[ir]
        inputTJLF.Q_LOC = abs(inputsPR.Q[ir])
        inputTJLF.Q_PRIME_LOC = inputsPR.Q_PRIME[ir]

        sum0 = 0
        for i = 1:ns
            sum0 = sum0 + inputTJLF.AS[i]*inputTJLF.TAUS[i]*(inputTJLF.RLNS[i]+inputTJLF.RLTS[i])
        end
        sum1 = 0
        for i = 1:ns
            sum1 = sum1 + inputsPR.AS[ir, i]*inputsPR.TAUS[ir, i]*(inputsPR.RLNS[ir, i]+inputsPR.RLTS[ir, i])
        end
        inputTJLF.P_PRIME_LOC = inputsPR.P_PRIME[ir]*sum0/sum1
        
    end

    if (inputsPR.ROTATION_FLAG == 1)
        for i = 1:ns
            inputTJLF.VPAR[i] = inputsPR.VPAR[ir, i]
            inputTJLF.VPAR_SHEAR[i] = inputsPR.VPAR_SHEAR[ir, i]
        end
    else
        for i = 1:ns
            inputTJLF.VPAR[i] = 0.0
            inputTJLF.VPAR_SHEAR[i] = 0.0
        end
    end

    inputTJLF.USE_BPER = true
    inputTJLF.USE_BPAR = false

    inputTJLF.BETAE = inputsPR.BETAE[ir]
    inputTJLF.XNUE = 0.0
    inputTJLF.ZEFF = inputsPR.ZEFF[ir]


    if (inputsEP.MODE_IN == 4)
        inputTJLF.FILTER = 2.0
    end

    kym = inputsEP.KY_MODEL
    if (kym == 0)
        inputTJLF.KY = 0.01*inputsEP.NTOROIDAL
    elseif (kym == 1)
        inputTJLF.KY = inputsEP.NTOROIDAL*inputTJLF.Q_LOC/inputTJLF.RMIN_LOC*inputsPR.RHO_STAR[ir]
    elseif (kym == 2)
        inputTJLF.KY = inputsEP.NTOROIDAL*0.1*inputTJLF.ZS[is]/sqrt(inputTJLF.MASS[is]*inputTJLF.TAUS[is])
    elseif (kym == 3)
        # This depends on a previous definition in kwscale_scan...
        inputTJLF.KY = inputsEP.KYHAT_IN*inputTJLF.ZS[is]/sqrt(inputTJLF.MASS[is]*inputTJLF.TAUS[is])
    end

    # This is one of the only things that is ran to for inputTJLF:
    inputsEP.FREQ_AE_UPPER = -abs(TJLFEP.exproConst.omegaGAM[ir])
    if inputsEP.ROTATIONAL_SUPPRESSION_FLAG == 1
        inputsEP.GAMMA_THRESH_MAX = abs(TJLFEP.exproConst.gammap[ir]) * 2.0 * (min(1.0 - inputsPR.RMIN[ir], inputsPR.RMIN[ir]) / inputsPR.RMAJ[ir])
        inputsEP.GAMMA_THRESH = 0.15 * abs(TJLFEP.exproConst.gammaE[ir] / inputsPR.SHEAR[ir])   # Bass PoP 2017 flow-shear suppression of AEs
        inputsEP.GAMMA_THRESH = min(inputsEP.GAMMA_THRESH, inputsEP.GAMMA_THRESH_MAX)
    else
        
        inputsEP.GAMMA_THRESH = 1.0e-7
        inputsEP.GAMMA_THRESH_MAX = 1.0e-7
    end

    # println("gammaE", inputsPR.gammaE)
    # println("Gamma Thresh MAX", inputsEP.GAMMA_THRESH_MAX)
    for field in fieldnames(typeof(inputsEP))
        value = getfield(inputsEP, field)
        #println("$field: $value")
        # println(fieldnames(inputTJLF))
        # println(fieldnames(InputTJLF))
        # if value[1] ==NaN
        #     value = 0
        # elseif value[1]  ==missing
        #     value[1] == coalesce(input.TJLF, false)
    end
# Create an instance of your struct (replace with your actual struct's name)

    # for field in fieldnames(typeof(inputTJLF))

        # if inputTJLF.field[1]== NaN
        #     inputTJLF.field=tesinput.field
        
        # end

    # for field in fieldnames(typeof(inputTJLF))
    #     value = getfield(inputTJLF, field)
    #     println("$field: $value")
        
    # end
    return inputTJLF
end
"""
readEXPRO function is a temporary function just used to define any EXPRO constants that are needed.
It is not going to be used in driver or mainsub or such.
"""
function readEXPRO(filename::String, is_EP::Int64)
    #filename = "/Users/benagnew/TJLF.jl/outputs/tglfep_tests/input.EXPRO"
    open(filename)

    lines = readlines(filename)

    ni1::Vector{Float64} = fill(NaN, 201)
    ni2::Vector{Float64} = fill(NaN, 201)
    ni3::Vector{Float64} = fill(NaN, 201)
    ni4::Vector{Float64} = fill(NaN, 201)
    Ti1::Vector{Float64} = fill(NaN, 201)
    Ti2::Vector{Float64} = fill(NaN, 201)
    Ti3::Vector{Float64} = fill(NaN, 201)
    Ti4::Vector{Float64} = fill(NaN, 201)
    dlnnidr1::Vector{Float64} = fill(NaN, 201)
    dlnnidr2::Vector{Float64} = fill(NaN, 201)
    dlnnidr3::Vector{Float64} = fill(NaN, 201)
    dlnnidr4::Vector{Float64} = fill(NaN, 201)
    dlntidr1::Vector{Float64} = fill(NaN, 201)
    dlntidr2::Vector{Float64} = fill(NaN, 201)
    dlntidr3::Vector{Float64} = fill(NaN, 201)
    dlntidr4::Vector{Float64} = fill(NaN, 201)
    cs::Vector{Float64} = fill(NaN, 201)
    rmin_ex::Vector{Float64} = fill(NaN, 201)
    gammaE::Vector{Float64} = fill(NaN, 201)
    gammap::Vector{Float64} = fill(NaN, 201)
    omegaGAM::Vector{Float64} = fill(NaN, 201)

    for line in lines[1:length(lines)]
        line = split(line, "=")
        val = line[2]
        line = split(line[1], "_")
        for i in eachindex(line)
            line[i] = strip(line[i])
        end

        exproname = line[2]
        if (length(line) == 4)
            isEPname = line[3]
        elseif (length(line) == 3)
            isEPname = ""
        end
        index = parse(Int64, String(line[end]))
        name = exproname*isEPname
        if (name == "ni1")
            ni1[index] = parse(Float64, String(val))
        elseif (name == "ni2")
            ni2[index] = parse(Float64, String(val))
        elseif (name == "ni3")
            ni3[index] = parse(Float64, String(val))
        elseif (name == "ni4")
            ni4[index] = parse(Float64, String(val))
        elseif (name == "Ti1") 
            Ti1[index] = parse(Float64, String(val))
        elseif (name == "Ti2")
            Ti2[index] = parse(Float64, String(val))
        elseif (name == "Ti3")
            Ti3[index] = parse(Float64, String(val))
        elseif (name == "Ti4")
            Ti4[index] = parse(Float64, String(val))
        elseif (name == "dlnnidr1")
            dlnnidr1[index] = parse(Float64, String(val))
        elseif (name == "dlnnidr2")
            dlnnidr2[index] = parse(Float64, String(val))
        elseif (name == "dlnnidr3")
            dlnnidr3[index] = parse(Float64, String(val))
        elseif (name == "dlnnidr4")
            dlnnidr4[index] = parse(Float64, String(val))
        elseif (name == "dlntidr1")
            dlntidr1[index] = parse(Float64, String(val))
        elseif (name == "dlntidr2")
            dlntidr2[index] = parse(Float64, String(val))
        elseif (name == "dlntidr3")
            dlntidr3[index] = parse(Float64, String(val))
        elseif (name == "dlntidr4")
            dlntidr4[index] = parse(Float64, String(val))
        elseif (name == "cs")
            cs[index] = parse(Float64, String(val))
        elseif (name == "rmin")
            rmin_ex[index] = parse(Float64, String(val))
        elseif (name == "gammaE")
            gammaE[index] = parse(Float64, String(val))
        elseif (name == "gammap")
            gammap[index] = parse(Float64, String(val))
        elseif (name == "omegaGAM")
            omegaGAM[index] = parse(Float64, String(val))
        end
    end
    # println("gammaE", gammaE)
    # println("gammap", gammap)
    # println("omegaGAM", omegaGAM)
    # println("cs", cs)
    # println("rmin", rmin_ex)
    # Diverge 4 is_EP values for each quatnity:

    if (is_EP == 1)
        return ni1, Ti1, dlnnidr1, dlntidr1, cs, rmin_ex, gammaE, gammap, omegaGAM
    elseif (is_EP == 2)
        return ni2, Ti2, dlnnidr2, dlntidr2, cs, rmin_ex, gammaE, gammap, omegaGAM
    elseif (is_EP == 3)
        return ni3, Ti3, dlnnidr3, dlntidr3, cs, rmin_ex, gammaE, gammap, omegaGAM
    elseif (is_EP == 4)
        return ni4, Ti4, dlnnidr4, dlntidr4, cs, rmin_ex, gammaE, gammap, omegaGAM
    else
        println("is_EP not within range. Check input.TGLFEP input")
        return 1
    end

end 
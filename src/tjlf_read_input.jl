"""
    readInput(baseDirectory::String)

parameter:
    baseDirectory::String - string of the directory (include a '\\' at the end) of input.tglf (NOT input.tglf.gen)

return:
    inputTJLF::InputTJLF - return InputTJLF struct based off the input.tglf file

description:
    parse through a input.tglf file found in baseDirectory parameter, creates a inputTJLF struct and populates the fields
    based off the values in the file. has some check to make sure the file is written properly
"""
function readInput(filename::String)::InputTJLF
    # gets the input.tglf file
    lines = readlines(filename)

    # finds the # of species to create Vectors later
    ns = -1
    nky = -1
    kygrid_model = -1
    for line in lines[1:length(lines)]
        line = split(line, "\n")
        line = split(line[1],"=")
        if line[1] == "NS"
            ns = parse(Int, strip(line[2]))
        elseif line[1] == "NKY"
            nky = parse(Int, strip(line[2]))
        elseif line[1] == "KYGRID_MODEL"
            kygrid_model = parse(Int, strip(line[2]))
        end
    end
    # make sure ns is defined
    @assert ns!=-1 "did not find NS in $filename (make sure this is an input.tglf file)"
    @assert nky!=-1 "did not find NKY"
    @assert kygrid_model!=-1 "did not find KYGRID_MODEL"
    @assert kygrid_model>=0 && kygrid_model <=5 "KYGRID_MODEL must be Int between 0 and 5"
    nky = get_ky_spectrum_size(nky,kygrid_model)



    # create InputTJLF struct
    inputTJLF = InputTJLF{Float64}(ns,nky)
    # fields that aren't used or implemented
    deletedFields = ["USE_TRANSPORT_MODEL","GEOMETRY_FLAG","B_MODEL_SA","FT_MODEL_SA","VPAR_SHEAR_MODEL","WRITE_WAVEFUNCTION_FLAG","VTS_SHEAR","VNS_SHEAR","VEXB","RMIN_SA","RMAJ_SA","Q_SA","SHAT_SA","ALPHA_SA","XWELL_SA","THETA0_SA","NN_MAX_ERROR"]

    # go through each line of the input.tglf file
    for line in lines[1:length(lines)]
        if contains(line,'#') || !contains(line,'=') continue end # ignores comments in the input.tglf file
        # gets rid of white space
        line = replace(line, " "=>"")
        # splits the line so that line[1] is the variable name and line[2] is the variable value
        line = split(line, "\n")
        line = split(line[1],"=")

        if line[1] ∈ deletedFields continue end 

        # check for the species field since they are all named [field name]_[species number]
        check = match(r"_\d",line[1])
        
        if check !== nothing # is a species field

            if replace(line[1], r"_\d"=>"") ∈ deletedFields continue end
            # alternative way of getting field nad species number
            # temp = split(line[1],"_")
            # get field and species number
            speciesField = Symbol(replace(line[1], r"_\d"=>""))
            speciesIndex = check.match[2:end]
            # skip values beyond the species number
            if parse(Int,speciesIndex) > ns continue end
            # set the value to the species field vectors
            getfield(inputTJLF, speciesField)[parse(Int,speciesIndex)] = parse(Float64,strip(line[2], ['\'','.',' ']))

        # species vector as a vector
        elseif line[2][1] == '['
            field = Symbol(line[1])
            getfield(inputTJLF, field) .= [parse(Float64,item) for item in split(line[2][2:end-1],",")]
        
        else # if not for the species vector

            # get the field name
            field = Symbol(line[1])

            # string
            if line[2][1] == '\'' || line[2][1] == '\"'
                val = string(strip(line[2], ['\'']))

            # bool
            elseif line[2][1] == '.'
                # make sure the input file is written correctly
                if lowercase(strip(line[2], ['\'','.'])) != "true" && lowercase(strip(line[2], ['\'','.'])) != "false"
                    println(line[1])
                    error("please put quotes around your string")
                end
                val = lowercase(strip(line[2], ['\'','.'])) == "true"

            # int
            elseif !contains(line[2],'.')
                val = parse(Int, strip(line[2], ['\'','.',' ']))

            # float
            else
                val = parse(Float64,strip(line[2], ['\'','.',' ']))
            end

            # set the inputTJLF field value
            try
                setfield!(inputTJLF,field,val)
            catch
                println(field)
                println(val)
                throw(error(field))
            end

        end
    end

    # just hard coded for some reason lol
    if inputTJLF.SAT_RULE == 2 || inputTJLF.SAT_RULE == 3
        inputTJLF.UNITS = "CGYRO"
        ### WTF
        inputTJLF.XNU_MODEL = 3
        inputTJLF.WDIA_TRAPPED = 1.0
    end

    inputTJLF.WIDTH_SPECTRUM .= inputTJLF.WIDTH
    inputTJLF.KY_SPECTRUM .= NaN
    inputTJLF.EIGEN_SPECTRUM .= NaN

    # double check struct is properly populated

    #If you want to test a long-format example (see tjlf_modules.jl), inP = true, otherwise, inP = false.
    #This is overly simplistic but I didn't want to mess anything else up.
    inP = false
    if inP
        inputTJLF = InputTJLF{Float64}(inP)
    end
    
    #Maybe checkInput could be altered for inputting default values, or the struct in modules could be
    #Redefined as having 
    checkInput(inputTJLF)

    return inputTJLF

end


"""
    checkInput(inputTJLF::InputTJLF)

description:
    check that the InputTJLF struct is properly populated
"""
function checkInput(inputTJLF::InputTJLF)
    field_names = fieldnames(InputTJLF)
    for field_name in field_names
        field_value = getfield(inputTJLF, field_name)
        if typeof(field_value)<:Missing
            @assert !ismissing(field_value) "Did not properly populate inputTJLF for $field_name = $field_value"
        end
        if typeof(field_value)<:Real
            @assert !isnan(field_value) "Did not properly populate inputTJLF for $field_name = $field_value"
        end
        if typeof(field_value)<:Vector && field_name!=:KY_SPECTRUM && field_name!=:EIGEN_SPECTRUM
            for val in field_value
                @assert !isnan(val) "Did not properly populate inputTJLF for array $field_name = $val"
            end
        end
    end
    if !inputTJLF.FIND_EIGEN
        @assert !inputTJLF.FIND_WIDTH "If FIND_EIGEN false, FIND_WIDTH should also be false"
    end
end

function checkInput(inputTJLFVector::Vector{InputTJLF})
    for inputTJLF in inputTJLFVector
        field_names = fieldnames(InputTJLF)
        for field_name in field_names
            field_value = getfield(inputTJLF, field_name)
            if typeof(field_value)<:Missing
                @assert !ismissing(field_value) "Did not properly populate inputTJLF for $field_name = $field_value"
            end
            if typeof(field_value)<:Real
                @assert !isnan(field_value) "Did not properly populate inputTJLF for $field_name = $field_value"
            end
            if typeof(field_value)<:Vector && field_name!=:KY_SPECTRUM && field_name!=:EIGEN_SPECTRUM
                for val in field_value
                    @assert !isnan(val) "Did not properly populate inputTJLF for array $field_name = $val"
                end
            end
        end
        if !inputTJLF.FIND_EIGEN
            @assert !inputTJLF.FIND_WIDTH "If FIND_EIGEN false, FIND_WIDTH should also be false"
        end
    end
end
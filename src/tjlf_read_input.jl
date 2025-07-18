"""
    readInput(baseDirectory::String)

parse through a input.tglf file found in baseDirectory parameter, creates a inputTJLF struct and populates the fields
based off the values in the file. has some check to make sure the file is written properly

parameter:

  - baseDirectory::String - string of the directory (include a '\\' at the end) of input.tglf (NOT input.tglf.gen)

return:

  - inputTJLF::InputTJLF - return InputTJLF struct based off the input.tglf file
"""
function readInput(filename::String)::InputTJLF
    # gets the input.tglf file
    lines = readlines(filename)

    # finds the # of species to create Vectors later
    ns = -1
    nky = -1
    kygrid_model = -1
    for line in lines
        if contains(line, '#') || !contains(line, '=')
            continue
        end # ignores comments in the input.tglf file
        field, value = strip.(split(line, "="))
        if field == "NS"
            ns = parse(Int, strip(value))
        elseif field == "NKY"
            nky = parse(Int, strip(value))
        elseif field == "KYGRID_MODEL"
            kygrid_model = parse(Int, strip(value))
        end
    end
    # make sure ns is defined
    @assert ns != -1 "did not find NS in $filename (make sure this is an input.tglf file)"
    @assert nky != -1 "did not find NKY"
    @assert kygrid_model != -1 "did not find KYGRID_MODEL"
    @assert kygrid_model >= 0 && kygrid_model <= 5 "KYGRID_MODEL must be Int between 0 and 5"
    nky = get_ky_spectrum_size(nky, kygrid_model)

    # create InputTJLF struct
    inputTJLF = InputTJLF{Float64}(ns, nky)
    # fields that aren't used or implemented
    deletedFields = [
        "USE_TRANSPORT_MODEL",
        "GEOMETRY_FLAG",
        "B_MODEL_SA",
        "FT_MODEL_SA",
        "VPAR_SHEAR_MODEL",
        "WRITE_WAVEFUNCTION_FLAG",
        "VTS_SHEAR",
        "VNS_SHEAR",
        "VEXB",
        "RMIN_SA",
        "RMAJ_SA",
        "Q_SA",
        "SHAT_SA",
        "ALPHA_SA",
        "XWELL_SA",
        "THETA0_SA",
        "NN_MAX_ERROR"
    ]

    # go through each line of the input.tglf file
    for line in lines
        if contains(line, '#') || !contains(line, '=')
            continue
        end # ignores comments in the input.tglf file

        # splits the line so that field is the variable name and value is the variable value
        field, value = strip.(split(line, "="))

        if field ∈ deletedFields
            continue
        end

        # check for the species field since they are all named [field name]_[species number]
        if match(r"_\d", field) !== nothing # is a species field
            if replace(field, r"_\d" => "") ∈ deletedFields
                continue
            end

            # get field and species number
            speciesField, speciesIndex = rsplit(field, "_"; limit=2)
            speciesField = Symbol(speciesField)
            speciesIndex = parse(Int, speciesIndex)
            # skip values beyond the species number
            if speciesIndex > ns
                continue
            end
            # set the value to the species field vectors
            getfield(inputTJLF, speciesField)[speciesIndex] = parse(Float64, value)

            # species vector as a vector
        elseif startswith(value, '[') 
            try
              getfield(inputTJLF, Symbol(field)) .= [parse(Float64, strip(item)) for item in split(value[2:end-1], ",")]
            catch e
              continue   # to skip NaN vectors KY_SPECTRUM
            end
        elseif startswith(value, "ComplexF64[")
            try 
             getfield(inputTJLF, Symbol(field)) .= [parse(ComplexF64, strip(item)) for item in split(value[12:end-1], ",")]
            catch e
              continue
            end
        else # if not for the species vector
            # string (quoted)
            if startswith(value, '\'') || startswith(value, '\"')
                val = string(strip(value, ['\'', '"']))

                # bool (handle both .true./.false. and true/false formats)
            elseif value == ".true." || value == ".false." || value == "true" || value == "false"
                val = (lowercase(value) == ".true." || lowercase(value) == "true")

                # int (only if it can actually be parsed as int)
            elseif tryparse(Int, value) !== nothing
                val = parse(Int, value)

                # float (only if it can actually be parsed as float)
            elseif tryparse(Float64, value) !== nothing
                val = parse(Float64, value)

                # string (fallback for anything else)
            else
                val = string(value)
            end

            # set the inputTJLF field value
            try
                setfield!(inputTJLF, Symbol(field), val)
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

    # Handle USE_BPER setting (consistent with TurbulentTransport.apply_presets!)
    if inputTJLF.USE_BPER
        inputTJLF.ALPHA_MACH = 0.0
    end

    inputTJLF.WIDTH_SPECTRUM .= inputTJLF.WIDTH
    if ismissing(inputTJLF.KY_SPECTRUM)
        inputTJLF.KY_SPECTRUM .= NaN 
    end
    if ismissing(inputTJLF.EIGEN_SPECTRUM)
        inputTJLF.EIGEN_SPECTRUM .= NaN 
    end
    

    # double check struct is properly populated

    if ismissing(inputTJLF.FIND_EIGEN)
        inputTJLF.FIND_EIGEN = true
    end

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

check that the InputTJLF struct is properly populated
"""
function checkInput(inputTJLF::InputTJLF)
    field_names = fieldnames(InputTJLF)
    for field_name in field_names
        field_value = getfield(inputTJLF, field_name)
        if typeof(field_value) <: Missing
            @assert !ismissing(field_value) "Did not properly populate inputTJLF for $field_name = $field_value"
        end
        if typeof(field_value) <: Real
            @assert !isnan(field_value) "Did not properly populate inputTJLF for $field_name = $field_value"
        end
        if typeof(field_value) <: Vector && field_name != :KY_SPECTRUM && field_name != :EIGEN_SPECTRUM
            for val in field_value
                @assert !isnan(val) "Did not properly populate inputTJLF for array $field_name = $field_value"
            end
        end
    end
    if !inputTJLF.FIND_EIGEN
        @assert !inputTJLF.FIND_WIDTH "If FIND_EIGEN false, FIND_WIDTH should also be false"
    end
end

function checkInput(inputTJLFVector::Vector{InputTJLF{T}}) where {T<:Real}
    for inputTJLF in inputTJLFVector
        checkInput(inputTJLF)
    end
end

"""
    save(input::InputTJLF, filename::AbstractString)

Write input_tjlf to file in InputTJLF format to be read by TJLF
"""
function save(input::InputTJLF, filename::AbstractString)
        open(filename, "w") do io
            for key in fieldnames(typeof(input))
                if startswith(String(key), "_")
                    continue
                end
                try
                    value = getfield(input, key)
                    if ismissing(value)
                        continue
                    elseif isa(value, Int)
                        println(io, "$(key)=$(convert(Int, value))")
                    elseif isa(value, String)
                        println(io, "$(key)='$value'")
                    elseif isa(value, Bool)
                        println(io, "$(key)=.$value.")
                    elseif isa(value, Vector{Float64})
                        println(io, "$(key)=$(convert(Vector{Float64}, value))")
                    elseif isa(value, Vector{ComplexF64})
                        println(io, "$(key)=$(convert(Vector{ComplexF64}, value))")
                    else
                        println(io, "$(key)=$(convert(Float64, value))")
                    end
                catch e
                    println("Error writing $key to input file")
                    rethrow(e)
                end
            end
        end
   
end


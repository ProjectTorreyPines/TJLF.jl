"""
    function readInput(baseDirectory::String)

parameter:
    baseDirectory::String - string of the directory (include a '\\' at the end) of input.tglf (NOT input.tglf.gen)

return:
    inputTJLF::InputTJLF - return InputTJLF struct based off the input.tglf file

description:
    parse through a input.tglf file found in baseDirectory parameter, creates a inputTJLF struct and populates the fields
    based off the values in the file. has some check to make sure the file is written properly

"""
function readInput(baseDirectory::String)::InputTJLF
    # gets the input.tglf file
    fileDirectory = joinpath(baseDirectory, "input.tglf")
    lines = readlines(fileDirectory)

    # finds the # of species to create Vectors later
    ns = -1
    nky = -1
    kygrid_model = -1
    nwidth = -1
    for line in lines[1:length(lines)]
        line = split(line, "\n")
        line = split(line[1],"=")
        if line[1] == "NS"
            ns = parse(Int, strip(line[2]))
        elseif line[1] == "NKY"
            nky = parse(Int, strip(line[2]))
        elseif line[1] == "KYGRID_MODEL"
            kygrid_model = parse(Int, strip(line[2]))
        elseif line[1] == "NWIDTH"
            nwidth = parse(Int, strip(line[2]))
        end
    end
    # make sure ns is defined
    @assert ns!=-1 "did not find NS"
    @assert nky!=-1 "did not find NKY"
    @assert kygrid_model!=-1 "did not find KYGRID_MODEL"
    @assert kygrid_model>=0 && kygrid_model <=5 "KYGRID_MODEL must be Int between 0 and 5"

    nky = get_ky_spectrum_size(nky,kygrid_model)
    ############# nky should equal NWIDTH ##################
    @assert nky == nwidth "ky spectrum size should equal width spectrum size"

    # create InputTJLF struct
    inputTJLF = InputTJLF{Float64}(ns,nwidth)

    # go through each line of the input.tglf file
    for line in lines[1:length(lines)]
        if contains(line,'#') || !contains(line,'=') continue end # ignores comments in the input.tglf file
        # gets rid of white space
        line = replace(line, " "=>"")
        # splits the line so that line[1] is the variable name and line[2] is the variable value
        line = split(line, "\n")
        line = split(line[1],"=")

        # check for the species field since they are all named [field name]_[species number]
        check = match(r"_\d",line[1])
        
        if check !== nothing # is a species field

            # alternative way of getting field nad species number
            # temp = split(line[1],"_")
            # get field and species number
            speciesField = Symbol(replace(line[1], r"_\d"=>""))
            speciesIndex = check.match[2:end]
            # skip values beyond the species number
            if parse(Int,speciesIndex) > ns continue end
            # set the value to the species field vectors
            getfield(inputTJLF, speciesField)[parse(Int,speciesIndex)] = parse(Float64,strip(line[2], ['\'','.',' ']))
        
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

    # double check struct is properly populated
    field_names = fieldnames(InputTJLF)
    for field_name in field_names
        field_value = getfield(inputTJLF, field_name)
        if typeof(field_value)<:Real
            @assert !isnan(field_value) && !ismissing(field_value) "Did not properly populate inputTJLF for $field_name"
        end
        if typeof(field_value)<:Vector && field_name!=:KY_SPECTRUM
            for val in field_value
                @assert !isnan(val) "Did not properly populate inputTJLF for array $field_name"
            end
        end
    end

    return inputTJLF

end
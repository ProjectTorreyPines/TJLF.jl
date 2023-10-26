function readInput(baseDirectory::String)
    fileDirectory = joinpath(baseDirectory, "input.tglf")
    lines = readlines(fileDirectory)

    ns = -1
    for line in lines[1:length(lines)]
        line = split(line, "\n")
        line = split(line[1],"=")
        if line[1] == "NS"
            ns = parse(Int, strip(line[2]))
            break
        end
    end
    @assert ns!=-1

    inputTJLF = InputTJLF{Float64}(ns)

    for line in lines[1:length(lines)]
        if contains(line,'#') || !contains(line,'=') continue end
        line = replace(line, " "=>"")
        line = split(line, "\n")
        line = split(line[1],"=")

        #### for the species vector
        check = match(r"_\d",line[1])
        if check !== nothing
            temp = split(line[1],"_")
            speciesField = Symbol(replace(line[1], r"_\d"=>""))
            speciesIndex = check.match[2:end]
            if parse(Int,speciesIndex) > length(inputTJLF.ZS) continue end
            getfield(inputTJLF, speciesField)[parse(Int,speciesIndex)] = parse(Float64,strip(line[2], ['\'','.',' ']))
            # setfield!(inputSpecies[parse(Int,speciesIndex)],    speciesField,     parse(Float64,strip(line[2], ['\'','.',' '])))
        else # if not for the species vector
            field = Symbol(line[1])

            # string
            if line[2][1] == '\''
                val = string(strip(line[2], ['\'']))
            # bool
            elseif line[2][1] == '.'
                val = strip(line[2], ['\'','.']) == "true"
            # int
            elseif !contains(line[2],'.')
                val = parse(Int, strip(line[2], ['\'','.',' ']))
            # float
            else
                val = parse(Float64,strip(line[2], ['\'','.',' ']))
            end

            try
                setfield!(inputTJLF,field,val)
            catch
                println(field)
                println(val)
                throw(error(field))
            end

        end
    end

    if inputTJLF.SAT_RULE == 2 || inputTJLF.SAT_RULE == 3
        inputTJLF.UNITS = "CGYRO"
        ####### WTF
        inputTJLF.XNU_MODEL = 3
        inputTJLF.WDIA_TRAPPED = 1.0
    end


    return inputTJLF

end
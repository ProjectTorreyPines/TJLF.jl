"""
    readInput(baseDirectory::String)

parse through a input.tglf file found in baseDirectory parameter, creates a inputTJLF struct and populates the fields
based off the values in the file. has some check to make sure the file is written properly

parameter:

  - baseDirectory::String - string of the directory (include a '\\' at the end) of input.tglf (NOT input.tglf.gen)

return:

  - inputTJLF::InputTJLF - return InputTJLF struct based off the input.tglf file
"""
# Unknown keys are a hard error, never a fallback: Fortran TGLF silently ignores
# unrecognized input.tglf content and runs its default (SAT0) template — a typo'd
# key then produces a plausible-looking wrong answer. TJLF refuses instead.
_reject_unknown_key(field, filename) = throw(ArgumentError(
    "unknown key in $filename: '$field' — TJLF rejects unrecognized inputs " *
    "(Fortran TGLF would silently ignore it and run with defaults)"))

function readInput(filename::String)::InputTJLF
    # gets the input.tglf file
    lines = readlines(filename)

    # finds the # of species to create Vectors later
    ns = -1
    nky = -1
    kygrid_model = -1
    sat_rule_found = false
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
        elseif field == "SAT_RULE"
            sat_rule_found = true
        end
    end
    # make sure ns is defined
    @assert ns != -1 "did not find NS in $filename (make sure this is an input.tglf file)"
    @assert nky != -1 "did not find NKY"
    @assert kygrid_model != -1 "did not find KYGRID_MODEL"
    @assert kygrid_model >= 0 && kygrid_model <= 5 "KYGRID_MODEL must be Int between 0 and 5"
    # SAT_RULE must be explicit: the Fortran fallback (unset -> SAT0) has silently
    # produced wrong-saturation-rule runs too many times to honor it here.
    @assert sat_rule_found "did not find SAT_RULE in $filename — TJLF requires it explicitly (Fortran silently defaults to SAT_RULE=0)"
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
            hasfield(InputTJLF{Float64}, speciesField) || _reject_unknown_key(field, filename)
            # skip values beyond the species number
            if speciesIndex > ns
                continue
            end
            # set the value to the species field vectors
            getfield(inputTJLF, speciesField)[speciesIndex] = parse(Float64, value)

            # species vector as a vector
        elseif startswith(value, '[') 
            hasfield(InputTJLF{Float64}, Symbol(field)) || _reject_unknown_key(field, filename)
            try
              getfield(inputTJLF, Symbol(field)) .= [parse(Float64, strip(item)) for item in split(value[2:end-1], ",")]
            catch e
              continue   # to skip NaN vectors KY_SPECTRUM
            end
        elseif startswith(value, "ComplexF64[")
            hasfield(InputTJLF{Float64}, Symbol(field)) || _reject_unknown_key(field, filename)
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
            hasfield(InputTJLF{Float64}, Symbol(field)) || _reject_unknown_key(field, filename)
            try
                setfield!(inputTJLF, Symbol(field), val)
            catch
                throw(ArgumentError("bad value for $field in $filename: '$value' parsed as " *
                                    "$(typeof(val)), field expects $(fieldtype(InputTJLF{Float64}, Symbol(field)))"))
            end

        end
    end

    # validate the parsed UNITS before apply_presets! can overwrite it: SAT0 presets
    # force GYRO, which would otherwise silently paper over a typo'd UNITS line
    @assert inputTJLF.UNITS in ("GYRO", "CGYRO") "UNITS must be \"GYRO\" or \"CGYRO\" (got \"$(inputTJLF.UNITS)\" in $filename)"

    apply_presets!(inputTJLF)

    inputTJLF.WIDTH_SPECTRUM .= inputTJLF.WIDTH
    # KY_SPECTRUM/EIGEN_SPECTRUM are NaN-filled by the (ns,nky) constructor; if a caller
    # constructed via @kwdef defaults they are empty, in which case nothing to fill.
    if isempty(inputTJLF.KY_SPECTRUM)
        inputTJLF.KY_SPECTRUM = fill(NaN, length(inputTJLF.WIDTH_SPECTRUM))
    end
    if isempty(inputTJLF.EIGEN_SPECTRUM)
        inputTJLF.EIGEN_SPECTRUM = fill(ComplexF64(NaN), length(inputTJLF.WIDTH_SPECTRUM))
    end

    # double check struct is properly populated (FIND_EIGEN is a concrete Bool, default true)

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
    # value-domain checks: ints/bools/strings have no NaN sentinel, so catch the
    # common "unset or typo'd switch" cases here instead of silently misbehaving
    @assert inputTJLF.SAT_RULE in (0, 1, 2, 3) "SAT_RULE must be 0, 1, 2, or 3 (got $(inputTJLF.SAT_RULE))"
    @assert inputTJLF.UNITS in ("GYRO", "CGYRO") "UNITS must be \"GYRO\" or \"CGYRO\" (got \"$(inputTJLF.UNITS)\")"
    # only reachable with USE_PRESETS=false — SAT2/3 are defined in CGYRO units only
    if inputTJLF.SAT_RULE in (2, 3)
        @assert inputTJLF.UNITS == "CGYRO" "SAT_RULE=$(inputTJLF.SAT_RULE) requires UNITS=\"CGYRO\" (got \"$(inputTJLF.UNITS)\" with USE_PRESETS=false)"
    end
end

"""
    apply_presets!(inputTJLF::InputTJLF) -> InputTJLF

SAT_RULE-calibrated switch coupling, mirroring Fortran `tglf_startup.f90` (where the
gating `USE_PRESETS` is hard-coded `.TRUE.`; TJLF exposes it as an input field so the
coupling can be disabled deliberately). Runs on `readInput` and before every solve:

- SAT_RULE 2/3: `XNU_MODEL=3`, `WDIA_TRAPPED=1.0`, `UNITS` GYRO -> CGYRO
- SAT_RULE 1:   `XNU_MODEL=2`, `WDIA_TRAPPED=0.0`
- SAT_RULE 0:   `XNU_MODEL=2`, `WDIA_TRAPPED=0.0`, `UNITS="GYRO"`
- `USE_BPER=true`: `ALPHA_MACH=0.0`

These are part of each saturation rule's original calibration, so overriding them
requires `USE_PRESETS=false` (in which case checkInput still rejects SAT2/3+GYRO).
"""
function apply_presets!(inputTJLF::InputTJLF)
    inputTJLF.USE_PRESETS || return inputTJLF
    # Fortran resets wdia_trapped unconditionally under presets, then SAT2/3 raises it
    inputTJLF.WDIA_TRAPPED = 0.0
    if inputTJLF.SAT_RULE == 2 || inputTJLF.SAT_RULE == 3
        inputTJLF.XNU_MODEL = 3
        inputTJLF.WDIA_TRAPPED = 1.0
        if inputTJLF.UNITS == "GYRO"
            inputTJLF.UNITS = "CGYRO"
        end
    elseif inputTJLF.SAT_RULE == 1
        inputTJLF.XNU_MODEL = 2
    elseif inputTJLF.SAT_RULE == 0
        inputTJLF.UNITS = "GYRO"
        inputTJLF.XNU_MODEL = 2
    end
    if inputTJLF.USE_BPER
        inputTJLF.ALPHA_MACH = 0.0
    end
    return inputTJLF
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
                    if is_unset(value)   # NaN sentinel = never populated; don't write it
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


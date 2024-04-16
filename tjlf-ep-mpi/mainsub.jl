using MPI
function mainsub(inputsEP::InputTJLFEP, inputsPR::profile, COMM_IN::MPI.Comm, np_in::Int64, id_in::Int64)
    TJLFEP_COMM = COMM_IN

    id = id_in
    np = np_in
    x = inputsEP.PROCESS_IN
    #println("PROCESS_IN = ")
    #println(x)
    if (x == 1)
        msg = "No"
        return msg
    elseif (x == 2)
        msg = "No"
        return msg
    elseif (x == 3)
        msg = "No"
        return msg
    elseif (x == 4)
        # Will actually be done
        msg = "No"
        return msg
    elseif (x == 5)
        # Actively working on
        inputsEP.WIDTH_IN_FLAG = false
        inputsEP.MODE_IN = 2
        inputsEP.KY_MODEL = 3

        growthrate, inputsEP, inputsPR = kwscale_scan(inputsEP, inputsPR, TJLFEP_COMM, np, id)
        return growthrate, inputsEP, inputsPR
    elseif (x == 6)
        msg = "No"
        return msg # There will be a specific return statement for each process_in case.
    end


end
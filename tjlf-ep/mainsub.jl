using MPI
function mainsub(inputsEP::InputTJLFEP, inputsPR::profile, COMM_IN::MPI.Comm)
    TJLFEP_COMM = COMM_IN

    id = MPI.Comm_rank(TJLFEP_COMM)
    np = MPI.Comm_size(TJLFEP_COMM)
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

        growthrate = kwscale_scan(inputsEP, inputsPR, TJLFEP_COMM)
        return growthrate
        msg = "Yes"
        return msg
    elseif (x == 6)
        msg = "No"
        return msg # There will be a specific return statement for each process_in case.
    end


end
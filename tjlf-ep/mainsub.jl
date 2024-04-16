function mainsub(inputsEP::InputTJLFEP, inputsPR::profile)
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

        println(inputsEP.MODE_IN)

        growthrate, inputsEP, inputsPR = kwscale_scan(inputsEP, inputsPR)
        return growthrate, inputsEP, inputsPR
    elseif (x == 6)
        msg = "No"
        return msg # There will be a specific return statement for each process_in case.
    end
end
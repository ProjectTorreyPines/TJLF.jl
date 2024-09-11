function mainsub(inputsEP::InputTJLFEP, inputsPR::profile, printout::Bool = true)
    x = inputsEP.PROCESS_IN
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
        msg = "No"
        return msg
    elseif (x == 5)
        inputsEP.WIDTH_IN_FLAG = false
        # inputsEP.MODE_IN = 2
        inputsEP.KY_MODEL = 3

        if (printout)
            println(inputsEP.MODE_IN)
        end

        growthrate, inputsEP, inputsPR = kwscale_scan(inputsEP, inputsPR, printout)
        return growthrate, inputsEP, inputsPR
    elseif (x == 6)
        msg = "No"
        return msg 
    end
end
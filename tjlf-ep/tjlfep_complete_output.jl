"""
tjlfep_complete_output_complete_output is meant for completing the output for all
radii tested. 


"""
function tjlfep_complete_output(profile_in::Vector{T}, inputsEP::InputTJLFEP{Float64}, inputsPR::profile{Float64}) where (T<:Real)

    l_accept_profile = fill(true, inputsEP.SCAN_N)
    profile_out = fill(0.0, inputsPR.NR)

    ir_out = fill(0, inputsEP.SCAN_N)

    for i_r = inputsEP.IRS:inputsEP.IRS+inputsEP.SCAN_N-1
        ir_in = Int(i_r-inputsEP.IRS+1)
        if (inputsEP.INPUT_PROFILE_METHOD == 2)
            ir_out[ir_in] = Int(inputsEP.IR_EXP[ir_in])
        else
            ir_out[ir_in] = i_r
        end
        l_exclude = (profile_in[ir_in] == profile_in[ir_in]+1.0) # Infinity test (check validity)
        l_exclude = (l_exclude || (profile_in[ir_in] != profile_in[ir_in])) # NaN test (check validity)
        l_exclude = (l_exclude || (profile_in[ir_in] < 0.0)) # Negative value test (check validity)
        if (!l_exclude)
            profile_out[ir_out[ir_in]] = profile_in[ir_in]
        else
            profile_out[ir_out[ir_in]] = 0.0
            l_accept_profile[ir_in] = false
        end
    end # irs <= i_r <= ir_max_0

    for i_r = ir_out[inputsEP.SCAN_N]+1:inputsPR.NR
        profile_out[i_r] = 0.0
    end

    profile_mask = fill(1.0, inputsPR.NR)
    #

    for i_r = 2:inputsEP.SCAN_N-1
        if (!l_accept_profile[i_r-1] && !l_accept_profile[i_r+1])
            l_accept_profile[i_r] = false
            profile_mask[ir_out[i_r]] = 0.0
        end
    end

    for i_r in eachindex(profile_out)
        profile_out[i_r] = profile_mask[i_r]*profile_out[i_r]
    end


    ir_mark = fill(0, inputsPR.NR)
    ir_counter = 0

    for i_r = 1:inputsPR.NR
        if (profile_out[i_r] != 0.0)
            ir_counter += 1
            ir_mark[ir_counter] = i_r
        end
    end

    ir_counter_max = max(1, ir_counter)

    ir_min = ir_mark[1]
    ir_max = ir_mark[ir_counter_max]

    for i = 1:ir_mark[1]
        profile_out[i] = profile_out[ir_mark[1]]
    end

    ir_counter = 1

    for i_r = ir_mark[1]+1:ir_mark[ir_counter_max]
        ir_1 = ir_mark[ir_counter]
        ir_2 = ir_mark[ir_counter+1]
        if (profile_out[i_r] == 0.0)
            profile_out[i_r] = profile_out[ir_1] + (profile_out[ir_2]-profile_out[ir_1])*(inputsPR.RMIN[i_r]-inputsPR.RMIN[ir_1])/(inputsPR.RMIN[ir_2]-inputsPR.RMIN[ir_1])
        else
            ir_counter += 1
        end
    end

    for i = ir_mark[ir_counter_max]:inputsPR.NR
        profile_out[i] = profile_out[ir_mark[ir_counter_max]]
    end

    return profile_in, profile_out, ir_min, ir_max, l_accept_profile
end
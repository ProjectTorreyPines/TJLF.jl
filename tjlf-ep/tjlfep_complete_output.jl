"""
tjlfep_complete_output

inputs: profile_in, inputsEP, inputsPR

outputs: profile_in, profile_out, ir_min, ir_max, l_accept_profile

tjlfep_complete_output performs the interpolation of a specified profile given a set of radial locations tested.
profile_in takes a vector of total ir values and outputs a 201 long vector with this information.
"""
function tjlfep_complete_output(profile_in::Vector{T}, inputsEP::Options{Float64}, inputsPR::profile{Float64}) where (T<:Real)

    # Default accept all scans and allocate profile_out and ir_out:
    l_accept_profile = fill(true, inputsEP.SCAN_N)
    profile_out = fill(0.0, inputsPR.NR)

    ir_out = fill(0, inputsEP.SCAN_N)

    # This loop goes from the starting radius in increments of 1 for however many scans there are:
    # 3 scans goes IRS, IRS+1, IRS+2 for example...
    for i_r = inputsEP.IRS:(inputsEP.IRS+inputsEP.SCAN_N-1)
        
        # ir_in goes from 1:SCAN_N
        ir_in = i_r-inputsEP.IRS+1
        
        # set ir_out to IR_EXP, essentially (or other method which I didn't use at all):
        if (inputsEP.INPUT_PROFILE_METHOD == 2)
            ir_out[ir_in] = inputsEP.IR_EXP[ir_in]
        else
            ir_out[ir_in] = i_r
        end

        # If the value of profile of concern is infinite, not a number, or negative, set it to be excluded in the profile:
        l_exclude = (profile_in[ir_in] == Inf) # Infinity test (check validity)
        l_exclude = (l_exclude || (profile_in[ir_in] == NaN)) # NaN test (check validity)
        l_exclude = (l_exclude || (profile_in[ir_in] < 0.0)) # Negative value test (check validity)

        # if we aren't excluding it, set the profile_out at the ir_exp value we're testing to its actual value,
        # otherwise, set it to 0 and allot l_accept_profile at that scan to false.
        if (!l_exclude)
            profile_out[ir_out[ir_in]] = profile_in[ir_in]
        else
            profile_out[ir_out[ir_in]] = 0.0
            l_accept_profile[ir_in] = false
        end
    end # irs <= i_r <= ir_max_0

    # For radii above the final scan_n's radius to 201, set the profile to 0
    # I think this is a relic of MPI:
    for i_r = ir_out[inputsEP.SCAN_N]+1:inputsPR.NR
        profile_out[i_r] = 0.0
    end

    # Set profile_mask for interpolation:
    profile_mask = fill(1.0, inputsPR.NR)

    # For intermediate (not edge) radii that were scanned,
    # if the profile is rejected for both of its neighbors,
    # also reject the one in between; mark the mask as 0.0 in these cases only at these radii.
    for i_r = 2:inputsEP.SCAN_N-1
        if (!l_accept_profile[i_r-1] && !l_accept_profile[i_r+1])
            l_accept_profile[i_r] = false
            profile_mask[ir_out[i_r]] = 0.0
        end
    end

    # For all 201 radii, multiply the mask by the profile_out.
    # This does mean that only the radii that were scanned and accepted will remain in profile_out
    # No interpolation has happened yet.
    for i_r in eachindex(profile_out)
        profile_out[i_r] = profile_mask[i_r]*profile_out[i_r]
    end

    
    ir_mark = fill(0, inputsPR.NR)
    ir_counter = 0

    # Find the locations of these accepted radii:
    for i_r = 1:inputsPR.NR
        if (profile_out[i_r] != 0.0)
            ir_counter += 1
            ir_mark[ir_counter] = i_r
        end
    end

    # Find the maximum index for imark.
    # If all radii were rejected, return 1.
    ir_counter_max = max(1, ir_counter)

    # Set min and max values of accepted radii; if none were accepted,
    # they will both be 0. If only one was accepted, they will both be
    # that radius.
    ir_min = ir_mark[1]
    ir_max = ir_mark[ir_counter_max]

    # For all radii until the first accepted radius, set the values
    # of the profile to be equal to that radius' profile value.
    for i = 1:ir_mark[1]
        profile_out[i] = profile_out[ir_mark[1]]
    end

    # Reset counter:
    ir_counter = 1

    # For radii after the initial one until the max, set first radius to first accepted
    # and second to neighbor above. Using this, if the value at the first radius above 
    # the minimum is 0 (default), set it to an interpolated value based on the values
    # of the nearest accepted radii. Otherwise, skip over the esablished values, and move to the next
    # range of accepted radii.
    for i_r = ir_mark[1]+1:ir_mark[ir_counter_max]
        ir_1 = ir_mark[ir_counter]
        ir_2 = ir_mark[ir_counter+1]
        if (profile_out[i_r] == 0.0)
            profile_out[i_r] = profile_out[ir_1] + (profile_out[ir_2]-profile_out[ir_1])*(inputsPR.RMIN[i_r]-inputsPR.RMIN[ir_1])/(inputsPR.RMIN[ir_2]-inputsPR.RMIN[ir_1])
        else
            ir_counter += 1
        end
    end

    # This loop above means that if more than one scan doesn't satisfy this condition in a 3 scan run, nothing will be interpolated.


    # This loop sets all values above the max accepted ir to that radii's profile value.
    # For one scan accepted, this sets all values above to its value and for no scans accepted,
    # it maintains all set to 0. For 3 scans, if 201 is accepted, nothing changes here.
    for i = ir_mark[ir_counter_max]:inputsPR.NR
        profile_out[i] = profile_out[ir_mark[ir_counter_max]]
    end

    # We return the profile_in which is unchanged, the interpolated and altered profile_out,
    # ir_min (radius of minimum accepted scan) and ir_max (radius of max accepted scan),
    # as well as the l_accept_profile for each scan.
    return profile_in, profile_out, ir_min, ir_max, l_accept_profile
end
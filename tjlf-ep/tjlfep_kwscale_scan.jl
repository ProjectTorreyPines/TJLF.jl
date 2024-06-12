function kwscale_scan(inputsEP::InputTJLFEP{Float64}, inputsPR::profile{Float64}, printout::Bool = true)
    # These are for testing purposes:
    #baseDirectory = "/Users/benagnew/TJLF.jl/outputs/tglfep_tests/input.MTGLF"
    #inputsPR = readMTGLF(baseDirectory)
    
    nfactor = 10
    nefwid = 10
    nkyhat = 5
    nkwf = nfactor*nefwid*nkyhat
    k_max = 4
    l_write_out = true

    #TJLFEP_COMM = COMM_IN
    #np = np_in
    # Here, I deal with the issue of np needing to be local in kwscale_scan:
    #np_local = MPI.Comm_size(TJLFEP_COMM)
    #id = id_in
    #println(np_local, " np_local at id & ir", id, " ", inputsEP.IR)

    # So in kwscale_scan, this should happen for each color/ir.
    # the k loop needs to happen sequentially, but the 1:500 inner one doesn't.
    # That one will be using a thread.

    kyhat_min = 0.0
    kyhat_max = 1.0

    growthrate = zeros(Float64, nkyhat, nefwid, nfactor, inputsEP.NMODES)
    growthrate_out = 0.0
    frequency = zeros(Float64, nkyhat, nefwid, nfactor, inputsEP.NMODES)
    frequency_out = 0.0
    
    lkeep_i = fill(true, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    ltearing_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_th_pinch_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_i_pinch_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_e_pinch_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_EP_pinch_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_max_outer_panel_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_QL_ratio_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))

    f0 = 0.0
    f1 = inputsEP.FACTOR_IN # 1 for first round
    w0 = inputsEP.WIDTH_MIN
    w1 = inputsEP.WIDTH_MAX
    kyhat_min = 0.0
    kyhat_max = 1.0
    kyhat0 = kyhat_min
    kyhat1 = kyhat_max
    factor = fill(NaN, nfactor)
    efwid = fill(NaN, nefwid)
    kyhat = fill(NaN, nkyhat)

    ikyhat_write = floor(Int, nkyhat/2) # 2
    #println(ikyhat_write)
    iefwid_write = floor(Int, nefwid/2) # 5
    ifactor_write = nfactor # 10

    # for scope purposes:
    inputTJLF = InputTJLF{Float64}(inputsPR.NS, 12, true)
    imark_min::Int64 = 0
    f_guess = fill(NaN, (nkyhat, nefwid))
    ikyhat_mark::Int64 = 0
    iefwid_mark::Int64 = 0
    for k = 1:k_max
        # k doesn't use Threads!
        factor .= NaN
        efwid .= NaN
        kyhat .= NaN
        for i = 1:nfactor
            factor[i] = ((f1-f0)/nfactor)*i+f0
            # k = 1: FACTOR_IN = 1.0
            # [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        end
        for i = 1:nefwid
            efwid[i] = ((w1-w0)/nefwid)*i+w0
        end
        for i = 1:nkyhat
            kyhat[i] = ((kyhat1-kyhat0)/nkyhat)*i+kyhat0
        end
        #MPI.Barrier(TJLFEP_COMM)
        #println(np_local)
        #=if (k == 2 && id == 0 && inputsEP.IR == 201 && false)
            println("id = 0 inputs before 2nd round:")
            println(inputsEP)
            #println(inputsPR)
        end
        #MPI.Barrier(TJLFEP_COMM)
        if (k == 2 && id == 1 && inputsEP.IR == 201 && false)
            println("id = 1 inputs before 2nd round:")
            println(inputsEP)
            #println(inputsPR)
        end
        if (inputsEP.IR == 201 && id == 0 && false)
            println("Round: ", k)
            println("(f0, f1), (id, ir): (", f0, ", ", f1, "), (", id, ", ", inputsEP.IR, ")")
            sleep(0.5)
            println("(w0, w1), (id, ir): (", w0, ", ", w1, "), (", id, ", ", inputsEP.IR, ")")
            sleep(0.5)
            println("(k0, k1), (id, ir): (", kyhat0, ", ", kyhat1, "), (", id, ", ", inputsEP.IR, ")")
            sleep(0.5)
            println("factor, efwid, kyhat: ", factor, efwid, kyhat)
        end=#
        # Threads:
        #=Threads.@threads=# 
        for i = 1:nkwf
            l_wavefunction_out = 0

            # The following 3 statments define each combination of ikyhat, iefwid, and ifactor.
            ikyhat = Int(floor((i-1)/(nefwid*nfactor))+1) # 2
            iefwid = Int(floor(1.0*mod(i-1, nefwid*nfactor)/nfactor)+1) # 5
            ifactor = mod(i-1, nfactor)+1 # 10
            if (k == 3 && i == 1)
                #println("=== efwid : iefwid : efwid[iefwid] ===")
                #println(efwid, " : ", iefwid, " : ", efwid[iefwid])
            end

            inputsEP.FACTOR_IN = factor[ifactor] # Just used in mapping
            inputsEP.KYHAT_IN = kyhat[ikyhat] # Set equal to ky_in
            inputsEP.WIDTH_IN = efwid[iefwid]
            #kyIndex = ikyhat
            #println("WIDTH_IN at pass ", i, ": ", inputsEP.WIDTH_IN)
            #println("np: ", np)

            #println(kyhat_in)
            #println(typeof(kyhat_in))

            #=if (id == 0 && inputsEP.IR == 201 && k == 1)
                println(np_local)
            end=#
            #MPI.Barrier(TJLFEP_COMM)
            #=if (id == 1 && inputsEP.IR == 201 && k == 1)
                println("WIDTH_IN for id = 1 in ir = 201: ", inputsEP.WIDTH_IN)
            end=#
            #=if (id == 0 && inputsEP.IR == 201 && k == 1)
                println("Before ky: ", inputsEP.L_TH_PINCH, " : ", i, " : ", k, " : ", id)
            end
            MPI.Barrier(TJLFEP_COMM)=#
            #=if (id == 1 && inputsEP.IR == 201 && k == 1 && i == 1)
                println("Before ky: ", inputsEP.L_TH_PINCH, " : ", i, " : ", k, " : ", id)
            end

                if ((id == 0 || id == 1) && inputsEP.IR == 201 && k == 2 && i == 1)
                    println("Before ky: ", inputsEP.L_TH_PINCH, " : ", i, " : ", k, " : ", id)
                    println(inputsEP)
                end=#
                #=MPI.Barrier(TJLFEP_COMM)
                if (id == 1 && inputsEP.IR == 201 && k == 1 && i == 1)
                    println(inputsEP.L_TH_PINCH, " : ", i, " : ", k)
                    println(inputsEP)
                end=#


            str_sf = string(Char(mod(floor(Int, inputsEP.FACTOR_IN/100.0), 10) + UInt32('0'))) *
                     string(Char(mod(floor(Int, inputsEP.FACTOR_IN/10.0), 10) + UInt32('0')))  *
                     string(Char(mod(floor(Int, inputsEP.FACTOR_IN), 10) + UInt32('0'))) *
                     "." *
                     string(Char(mod(floor(Int, 10*inputsEP.FACTOR_IN), 10) + UInt32('0'))) *
                     string(Char(mod(floor(Int, 100*inputsEP.FACTOR_IN), 10) + UInt32('0'))) *
                     string(Char(mod(floor(Int, 1000*inputsEP.FACTOR_IN), 10) + UInt32('0')))
            #=if (k == 2 && id == 0 && inputsEP.IR == 201 && i <= 5)
                #println("id = 0 str_sf: ")
                #println(str_sf)
                println("Input: ", inputsEP)
                println("iteration: ", i)
            end
            MPI.Barrier(TJLFEP_COMM)
            if (k == 2 && id == 1 && inputsEP.IR == 201 && i <= 5)
                #println("id = 1 str_sf: ")
                #println(str_sf)
                println("Input: ", inputsEP)
                println("iteration: ", i)
            end=#
            str_wf_file = "out.wavefunction"*inputsEP.SUFFIX*"_sf"*str_sf
            #=if (inputsEP.IR == 201 && k == 1)
                println("ikyhat_write for id in ir = 201: ", id, " ", ikyhat_write)
            end=#
            if ((inputsEP.WRITE_WAVEFUNCTION == 1) &&
                (ikyhat == ikyhat_write) &&
                (iefwid == iefwid_write) &&
                (ifactor == ifactor_write) &&
                (k == k_max) && printout)
                l_wavefunction_out = 1
                #println(str_sf, " for round ", k, ", id & ir: ", id, " ", inputsEP.IR)
            end
            
            if (inputsEP.IR == 3 && k == 1 && false)
                println("================= Iter: ", i, " ================")
                println(inputsEP.PROCESS_IN)
                println(inputsEP.THRESHOLD_FLAG)
                println(inputsEP.SCAN_METHOD)
                println(inputsEP.QL_THRESH_RATIO)
                println(inputsEP.Q_SCALE)
                println(inputsEP.KY_MODEL)
                println(inputsEP.FACTOR_IN)
                println(inputsEP.WIDTH_IN)
                println(inputsEP.KYHAT_IN)
                println(inputsEP.WIDTH_MIN)
                println(inputsEP.WIDTH_MAX)
                println(inputsEP.IS_EP)
                println(inputsEP.IR)
                println(inputsEP.JTSCALE)
                println(inputsEP.NN)
                println(inputsEP.FREQ_AE_UPPER)
            end

            #testid = false
            #if (inputsEP.IR == 101)
            #    testid = true
            #end

            if (inputsEP.IR == 101)
                #println("============== Iter: ", i)
            end

            gamma_out, freq_out, inputTJLF = TJLFEP_ky(inputsEP, inputsPR, str_wf_file, l_wavefunction_out, printout)
            
            #=if (id == 0 && inputsEP.IR == 201 && k == 1)
                println("After ky: ", inputsEP.L_TH_PINCH, " : ", i, " : ", k, " : ", id)
            end
            MPI.Barrier(TJLFEP_COMM)=#
                #=if (id == 1 && inputsEP.IR == 201 && k == 1 && i == 1)
                    println("After ky: ", inputsEP.L_TH_PINCH, " : ", i, " : ", k, " : ", id)
                end

                if ((id == 0 || id == 1) && inputsEP.IR == 201 && k == 2 && i == 1)
                    println("After ky: ", inputsEP.L_TH_PINCH, " : ", i, " : ", k, " : ", id)
                    println(inputsEP)
                end=#
            if (inputsEP.IR == 2 && k == 1)
                #println("gamma_out for iter on run 1: ", gamma_out, " ", i)
            end

            # it is performing the entire operation over all ky. I don't want to run this a bunch of times since it gets rid of the purposes
            # of using MPI. 
            #=if (l_wavefunction_out == 1)
                io10 = open("testoutput_"*str_wf_file, "w")
                println(io10, gamma_out)
                println(io10, freq_out)
                close(io10)
            end=#
            #=if (id == 0 && inputsEP.IR == 201 && k == 3)
                println("id 0 l_th_pinch at iter: ", inputsEP.L_TH_PINCH, ", ", i)
            end
            MPI.Barrier(TJLFEP_COMM)
            if (id == 1 && inputsEP.IR == 201 && k == 3)
                println("id 1 l_th_pinch at iter: ", inputsEP.L_TH_PINCH, ", ", i)
            end=#
            #Following this is the run of getting:
            if (inputsEP.IR == 101)
                #println("===========")
                #println("Before L_MAX_OUTER_PANEL: ")
                #println(inputsEP.L_MAX_OUTER_PANEL)
            end

            for n = 1:inputsEP.NMODES
                growthrate[ikyhat,iefwid,ifactor,n] = gamma_out[n]
                frequency[ikyhat,iefwid,ifactor,n]  = freq_out[n]
                lkeep_i[ikyhat,iefwid,ifactor,n] = inputsEP.LKEEP[n]
                ltearing_i[ikyhat,iefwid,ifactor,n] = inputsEP.LTEARING[n]
                l_th_pinch_i[ikyhat,iefwid,ifactor,n] = inputsEP.L_TH_PINCH[n]
                l_i_pinch_i[ikyhat,iefwid,ifactor,n] = inputsEP.L_I_PINCH[n]
                l_e_pinch_i[ikyhat,iefwid,ifactor,n] = inputsEP.L_E_PINCH[n]
                l_EP_pinch_i[ikyhat,iefwid,ifactor,n] = inputsEP.L_EP_PINCH[n]
                l_max_outer_panel_i[ikyhat,iefwid,ifactor,n] = inputsEP.L_MAX_OUTER_PANEL[n]
                l_QL_ratio_i[ikyhat,iefwid,ifactor,n] = inputsEP.L_QL_RATIO[n]
            end
            #println(growthrate)
            #println("Iteration ", i, ", id ", id)

            if (inputsEP.IR == 101 && printout)
                println("============== Iter: ", i)
                #println("growthrate, l_max_outer_panel_i at: [", ikyhat, ", ", iefwid, ", ", ifactor, "]")
                #println(growthrate[ikyhat, iefwid, ifactor, :])
                #println("-----After----")
                println(lkeep_i[ikyhat,iefwid,ifactor, :])
            end
            #sds = -1
            #@assert sds != -1 "End"
        end # end of MPI process collection
        
        # Ending here are this point for testing purposes:
        #if (id == 0 && inputsEP.IR == 201 && k == 1)
        #    println("GRate and ir = 201: ")
        #    println(growthrate)
        #end
        #MPI.Barrier(TJLFEP_COMM)
        #if (id == 0 && inputsEP.IR == 2 && k == 1)
        #    println("GRate and ir = 2: ")
        #    println(growthrate)
        #end
        # Allreduce performs the operation (2nd parameter) between each of the processes in a specified color Group
        # determined by TJLFEP_COMM, then distributes in back to each of the processes individually.
        #if (id == 0)
        #    println(growthrate[:, 1, 1, 1])
        #end
        #MPI.Barrier(TJLFEP_COMM)
        #MPI.Allreduce!(growthrate, MPI.SUM, TJLFEP_COMM)
        #MPI.Allreduce!(frequency, MPI.SUM, TJLFEP_COMM)
        #MPI.Allreduce!(lkeep_i, MPI.LAND, TJLFEP_COMM)
        #MPI.Allreduce!(ltearing_i, MPI.LOR, TJLFEP_COMM)
        #MPI.Allreduce!(l_th_pinch_i, MPI.LOR, TJLFEP_COMM)
        #MPI.Allreduce!(l_i_pinch_i, MPI.LOR, TJLFEP_COMM)
        #MPI.Allreduce!(l_e_pinch_i, MPI.LOR, TJLFEP_COMM)
        #MPI.Allreduce!(l_max_outer_panel_i, MPI.LOR, TJLFEP_COMM)
        #MPI.Allreduce!(l_QL_ratio_i, MPI.LOR, TJLFEP_COMM)
        #if (id == 0)
        #    println("growthrate for id = 0 and ir = ", inputsEP.IR, growthrate[1, 1, :, :])
        #end 
        #=if (id == 0 && inputsEP.IR == 201 && k == 1)
            println("lkeep_i for id 1: ")
            println(lkeep_i)
        end
        MPI.Barrier(TJLFEP_COMM)
        if (id == 1 && inputsEP.IR == 201 && k == 1)
            println("lkeep_i for id 0: ")
            println(lkeep_i)
        end=#
        # This loop creates a 5x10 matrix full of 11. It then runs
        # through all dimensions of lkeep_i, which is a reference matrix
        # telling you where each 
        imark = fill(nfactor+1, (nkyhat, nefwid))
        #if (id == 0)
        #    println(imark, " imark at id: ", id)
        #end

        # The following loop sets the matrix imark finds the first mode for a specific ikyhat, iefwid, ifactor combination
        # that is set to be kept in lkeep_i. This mode is thus representing that specific ifactor value for the run
        # and imark is set to ifactor for that ikyhat and iefwid. The same is also done for the factor values as once the first non-zeros
        # growthrate is found along the factor spectrum, we return back to the next width run. Essentially, we are trying to find
        # all of the first non-zero growthrates at each kyhat, width combo; these are then allotted to the imark matrix which
        # describes whether a kyhat, width combo even has a non-zero growthrate, and if it does, where in the scan did it find that.

        for ikyhat = 1:nkyhat
            for iefwid = 1:nefwid
                for ifactor = 1:nfactor
                    for n = 1:inputsEP.NMODES # For a specified mode, 
                        if (lkeep_i[ikyhat, iefwid, ifactor, n])
                            imark[ikyhat, iefwid] = ifactor
                            # these could range from 1 to 10.
                            break
                        end
                    end # n
                    if (imark[ikyhat,iefwid] <= nfactor) break end
                end # ifactor
            end # iefwid
        end #ikyhat
        #=if (k == 2 && id == 0 && inputsEP.IR == 201)
            println("imark for second round id = 0: ", imark)
        end
        sleep(1)
        if (k == 2 && id == 1 && inputsEP.IR == 201)
            println("imark for second round id = 1: ", imark)
        end
        #if (id == 0)
        #    println(imark, " imark at id: ", id)
        #end =#

        # This loop searches for the lowest value of imark for a specific
        # round of K.
        imark_min = nfactor + 1
        for ikyhat = 1:nkyhat
            for iefwid = 1:nefwid
                imark_min = min(imark[ikyhat, iefwid], imark_min)
            end
        end
        #if (id == 0)
        #    println(imark_min, " imark_min at id ", id)
        #end

        # Next we are going to set:
        # fmark - 
        # gmark - 
        # f_guess_mark -
        # gamma_mark_i_1 -
        # gamma_mark_i_2 -
        # f_mark_i - 
        # lkeep_ref - 
        
        fmark = 1.0e20
        gmark = 0.0
        f_guess_mark = 1.0e20
        gamma_mark_i_1 = fill(NaN, (nkyhat, nefwid))
        gamma_mark_i_2 = fill(NaN, (nkyhat, nefwid))
        f_mark_i = fill(NaN, (nkyhat, nefwid))
        lkeep_ref = fill(false, (nkyhat, nefwid))


        # if there are any unstable modes:
        if (imark_min <= nfactor)
            
            # default the marked values as the last ones
            ikyhat_mark = nkyhat
            iefwid_mark = nefwid

            # loop all combos of kyhat, width in imark:
            for ikyhat = 1:nkyhat
                for iefwid = 1:nefwid
                    imark_ref = nfactor-1 # ? What is the point ?
                    #println("pass 1: ", imark_ref)
                    
                    # if imark is not the final one of the factor spectrum, set the reference to imark
                    # and set it as kept; otherwise, if it's the last one, set the reference to imark-1
                    # (9) and set as kept, and if it's a default value of 11, set reference to 10.
                    if (imark[ikyhat, iefwid] < nfactor) # < 10
                        imark_ref = imark[ikyhat, iefwid]
                        #println("route 1: ", imark_ref)
                        lkeep_ref[ikyhat, iefwid] = true
                    else # 10 or 11
                        imark_ref = imark[ikyhat, iefwid] - 1
                        #println("route 2: ", imark_ref)
                        if (imark[ikyhat, iefwid] == nfactor) # == 10
                            lkeep_ref[ikyhat,iefwid] = true
                        end
                    end
                    # Now we set the f_g1 and f_g2 values, which allows us to guess the
                    # value of f_guess based on this equation. This equation is not clear yet to me.
                    f_g1 = factor[imark_ref]
                    # Set possible default values (as Fortran accesses it despite not being allotted):
                    if (imark_ref+1 < 11)
                        f_g2 = factor[imark_ref+1]
                    else
                        f_g2 = 0
                    end

                    # This is the portion that is NOT consistent with the Fortran and is
                    # causing problems due to default values:
                    gamma_g1 = -2.0
                    gamma_g2 = -1.0
                    
                    # For all modes, set gamma_g1 to the maximum growthrate of kept modes at this point of ikyhat, iefwid, and imark_ref
                    # This also sets gamma_g2 to the maximum neighboring factor value of growthrate (imark_ref+1).
                    for n = 1:inputsEP.NMODES
                        if (lkeep_i[ikyhat, iefwid, imark_ref, n])
                            gamma_g1 = max((growthrate[ikyhat, iefwid, imark_ref, n]), gamma_g1)
                        end
                        if (imark_ref+1 < 11)
                            if (lkeep_i[ikyhat, iefwid, imark_ref+1, n])
                                gamma_g2 = max((growthrate[ikyhat, iefwid, imark_ref+1, n]), gamma_g2)
                            end
                        end
                        if (imark_ref+1 == 11 && n == 4 && printout)
                            println("ikyhat, iefwid: ", ikyhat, " ", iefwid)
                            println("imarkref+1=11 : f_g1 & f_g2: ", f_g1, " ", f_g2)
                            println("imarkref+1=11 : gamma_g1 & gamma_g2: ", gamma_g1, " ", gamma_g2)
                        end
                    end
                    f_guess[ikyhat, iefwid] = (gamma_g1*f_g2-gamma_g2*f_g1)/(gamma_g1-gamma_g2)
                    gamma_mark_i_1[ikyhat, iefwid] = gamma_g1
                    gamma_mark_i_2[ikyhat, iefwid] = gamma_g2
                    f_mark_i[ikyhat, iefwid] = f_g1
                end # iefwid
            end # ikyhat
            # This next one seems to only happen on the lower of the TWO
            # ir for just the first round of k ???
            #println("=== lkeep_ref ===")
            #println(lkeep_ref)
            #println("Round: ", k)
            #println("=================")

            # Now we scan across imark again, looking for modes that are:
            # kept; have a max growthrate along the modes that is less than 95% of the neighboring maximum; have a factor value that is 
            # less than or equal to the current maximum... (starts at 1e20)
            # then if these are satisfied, find modes that:
            # have a growthrate that is larger than the current maximum (starts at 0) or
            # have a factor that is smaller than the max
            for ikyhat = 1:nkyhat
                for iefwid = 1:nefwid
                    if ((lkeep_ref[ikyhat, iefwid]) && (gamma_mark_i_1[ikyhat, iefwid] < 0.95*gamma_mark_i_2[ikyhat, iefwid]) &&
                        (f_mark_i[ikyhat, iefwid] <= fmark))
                        if ((gamma_mark_i_1[ikyhat, iefwid] > gmark) || (f_mark_i[ikyhat,iefwid] < fmark))
                            # If all of these are satisfied, set new marks away from default:
                            ikyhat_mark = ikyhat
                            ikyhat_write = ikyhat
                            iefwid_mark = iefwid
                            iefwid_write = iefwid
                            
                            if (inputsEP.IR == 101)
                                #println("Marks set: ", iefwid_mark, " : ", ikyhat_mark)
                                #println("conds: lkeep, gamma_mark_i_1, 0.95*gamma_mark_i_2, f_mark_i, gmark, fmark")
                                #println(lkeep_ref[ikyhat,iefwid], " ", gamma_mark_i_1[ikyhat, iefwid], " ", 0.95*gamma_mark_i_2[ikyhat,iefwid])
                                #println(f_mark_i[ikyhat,iefwid], " ", gmark, " ", fmark)
                            end

                            # then set new maximums (fmark and gmark)
                            # and set f_guess_mark which corresponds to a maximized growthrate and
                            # minimized factor.
                            gmark = gamma_mark_i_1[ikyhat, iefwid]
                            fmark = f_mark_i[ikyhat, iefwid]
                            f_guess_mark = f_guess[ikyhat, iefwid]

                            # This loops over all combos of kyhat and width to find a single guess for "f"
                            
                            if (inputsEP.IR == 101)
                                #println("update: ", gmark, " : ", fmark, " : ", f_guess_mark)
                                #println("===========")
                            end

                            #=if (id == 0)
                                println("Statements set! for id & ir and k: ", id, " ", inputsEP.IR, " ", k)
                                println("kywrite, wdwrite: ", ikyhat_mark, " ", iefwid_mark)
                                println("After: ", fmark, ", ", gmark, ", ", f_guess_mark)
                            end
                            sleep(0.25)
                            if (id == 1)
                                println("Statements set! for id & ir and k: ", id, " ", inputsEP.IR, " ", k)
                                println("kywrite, wdwrite: ", ikyhat_mark, " ", iefwid_mark)
                                println("After: ", fmark, ", ", gmark, ", ", f_guess_mark)
                            end=#
                        end
                    else
                        #println("Statements NOT set! for id & ir: ", id, " ", inputsEP.IR)
                    end
                end # iefwid
            end # ikyhat
        end # if ending (mode found)

        # Next is the writing of the out.scalefactor files.

        g = fill(NaN, inputsEP.NMODES)
        f = fill(NaN, inputsEP.NMODES)
        keep_label = fill("", inputsEP.NMODES)
        # These will be produced whenever I run from a specified directory, say
        # via a batch file or some execution command for TJLF-EP. Running by REPL
        # does create these files, but it places them in the TJLF.JL directory,
        # which isn't where I will want to store them in the end.
        if (l_write_out && printout)
            filename = "out.scalefactor"*inputsEP.SUFFIX
            iexist = isfile(filename)
            if (iexist)
                io = open(filename, "a")
            else
                io = open(filename, "w")
                println(io, "factor, (gamma(n), freq(n), flag, n=1, nmodes_in)")
                println(io, "flag key:  'K' mode is kept")
                if (inputsEP.REJECT_TEARING_FLAG == 1) println(io, "           'T' rejected for tearing parity")  end
                if (inputsEP.REJECT_I_PINCH_FLAG == 1) println(io, "           'Pi' rejected for ion pinch") end
                if (inputsEP.REJECT_E_PINCH_FLAG == 1) println(io, "           'Pe' rejected for electron pinch") end
                if (inputsEP.REJECT_TH_PINCH_FLAG == 1) println(io, "           'Pth' rejected for thermal pinch") end
                if (inputsEP.REJECT_EP_PINCH_FLAG == 1) println(io, "           'PEP' rejected for EP pinch") end
                if (inputsEP.REJECT_MAX_OUTER_PANEL_FLAG == 1) println(io, "           'OP' rejected for ballooning-space max outside first panel") end
                println(io, "           'QLR' rejected for QL ratio EP/|chi_i| < ", inputsEP.QL_THRESH_RATIO)
                println(io, "           'F' rejected for non-AE frequency > ", inputsEP.F_REAL[inputsEP.IR]*inputsEP.FREQ_AE_UPPER, " kHz")
                if (inputsEP.REAL_FREQ == 1) println(io, "Frequencies in real units, plasma frame [kHz]") end
            end
            
            ky_write = kyhat[ikyhat_write]*inputTJLF.ZS[inputsEP.IS_EP+1]/sqrt(inputTJLF.MASS[inputsEP.IS_EP+1]*inputTJLF.TAUS[inputsEP.IS_EP+1])

            println(io, "--------------- ky*rho_EP= ", kyhat[ikyhat_write], "  (ky*rho_s=", ky_write, ")", 
            "   width= ", efwid[iefwid_write], 
            "   scalefactor= ", fmark, " -------------------")

            for ifactor = 1:nfactor
                ikyhat = ikyhat_write
                iefwid = iefwid_write
                tlabelvec = fill("", inputsEP.NMODES)
                for n = 1:inputsEP.NMODES
                    g[n] = inputsEP.F_REAL[inputsEP.IR]*growthrate[ikyhat, iefwid, ifactor, n]
                    f[n] = inputsEP.F_REAL[inputsEP.IR]*frequency[ikyhat, iefwid, ifactor, n]
                    testlabel = ""
                    if (lkeep_i[ikyhat,iefwid,ifactor,n])
                        keep_label[n] = " K  "
                        testlabel = testlabel*" K "
                    elseif (ltearing_i[ikyhat,iefwid,ifactor,n])
                        keep_label[n] = " T  "
                        testlabel = testlabel*" T "
                    elseif (l_i_pinch_i[ikyhat,iefwid,ifactor,n])
                        keep_label[n] = " Pi "
                        testlabel = testlabel*" Pi "
                    elseif (l_e_pinch_i[ikyhat,iefwid,ifactor,n])
                        keep_label[n] = " Pe "
                        testlabel = testlabel*" Pe "
                    elseif (l_th_pinch_i[ikyhat,iefwid,ifactor,n])
                        keep_label[n] = " Pth"
                        testlabel = testlabel*" Pth "
                    elseif (l_EP_pinch_i[ikyhat,iefwid,ifactor,n])
                        keep_label[n] = " PEP"
                        testlabel = testlabel*" PEP "
                    elseif (l_max_outer_panel_i[ikyhat,iefwid,ifactor,n])
                        keep_label[n] = " OP"
                        testlabel = testlabel*" OP "
                    elseif (l_QL_ratio_i[ikyhat,iefwid,ifactor,n])
                        keep_label[n] = " QLR"
                        testlabel = testlabel*" QLR "
                    elseif (frequency[ikyhat,iefwid,ifactor,n] > inputsEP.FREQ_AE_UPPER)
                        keep_label[n] = " F  "
                        testlabel = testlabel*" F "
                    else
                        keep_label[n] = " ?  "
                        testlabel = testlabel*" ? "
                    end
                    tlabelvec[n] = testlabel
                end
                if (inputsEP.REAL_FREQ == 0)
                    for n = 1:inputsEP.NMODES
                        println(io, factor[ifactor], " ", g[n], " ", f[n], keep_label[n])
                        println(io, tlabelvec[n])
                    end
                else
                    for n = 1:inputsEP.NMODES
                        println(io, factor[ifactor], " ", g[n], " ", f[n], keep_label[n])
                        println(io, tlabelvec[n])
                    end
                end
            end # ifactor
            close(io)
        end 

        # After printing the info just calculated, if fmark is already significantly small, adjust 
        # ranges within the range, otherwise, move to a new order of 10 in the factor.
        if (fmark < 1.0e10)
            # If in the first scan round, set the new max factor to the minimized factor, otherwise, do
            # some more specific things to hone in each round (including in kyhat and width)
            if (k == 1)
                f1 = fmark
                f0 = 0.0
                #println("fmark > 1e10 : Pass ", k)
            else
                delf = (f1-f0)/3
                #=if (id == 0 && inputsEP.IR == 201)
                    println("delf, f1, f0 for id = 0: ", delf, ", ", f1, ", ", f0)
                end
                MPI.Barrier(TJLFEP_COMM)
                if (id == 1 && inputsEP.IR == 201)
                    println("delf, f1, f0 for id = 1: ", delf, ", ", f1, ", ", f0)
                end=#
                f1 = fmark + delf
                f0 = fmark - delf
                if (f0 < 0.0) 
                    f0 = 0.0 
                end
                delw = (w1-w0)/4.0
                w1 = efwid[iefwid_mark] + delw
                w0 = efwid[iefwid_mark] - delw
                if (w1 > inputsEP.WIDTH_MAX)
                    w0 = w0 - (w1-inputsEP.WIDTH_MAX)
                    w1 = inputsEP.WIDTH_MAX
                elseif (w0 < inputsEP.WIDTH_MIN)
                    w1 = w1 + (inputsEP.WIDTH_MIN-w0)
                    w0 = inputsEP.WIDTH_MIN
                end
                #println(iefwid_mark, " iefwid_mark")
                #println(efwid, " efwid")
                #println("out: (w1, w0) : delw = (", w1, ", ", w0, ") : ", delw)
                #println("fmark > 1e10 : Pass ", k)
                delky = (kyhat1-kyhat0)/4
                kyhat1 = kyhat[ikyhat_mark] + delky
                kyhat0 = kyhat[ikyhat_mark] - delky
                if (kyhat1 > kyhat_max)
                    kyhat0 = kyhat0 - (kyhat1-kyhat_max)
                    kyhat1 = kyhat_max
                elseif (kyhat0 < kyhat_min)
                    kyhat1 = kyhat1 + (kyhat_min-kyhat0)
                    kyhat0 = kyhat_min
                end
            end # k > 1
        else
            f0 = f1
            f1 = 10.0*f1
            #println("fmark < 1e10 : Pass ", k)
        end

    end # end of k
    #println(imark_min)
    
    # After the fourth round, if no unstable modes have been found, default the scalefactor which will be used In
    # the calculation in the driver to 10k, the width to its default (I think around 100) and same for kyhat.
    # Otherwise, set the factor to the best guess at ikyhat_mark, iefwid_mark and their corresponding width and kyhat
    # to those marks.
    if (imark_min > nfactor)
        # If, over the scan of k, there's no unstable modes, set each to lowest.
        inputsEP.FACTOR_IN = 10000
        inputsEP.WIDTH_IN = efwid[1]
        inputsEP.KYMARK = kyhat[1]

        if (printout)
            println("imark_min > nfactor ", inputsEP.FACTOR_IN, " ", inputsEP.WIDTH_IN, " ", inputsEP.KYMARK)
        end
    else
        # If, over the scan of k, there's an unstable mode, set each to each marked point.
        inputsEP.FACTOR_IN = f_guess[ikyhat_mark, iefwid_mark]
        inputsEP.WIDTH_IN = efwid[iefwid_mark]
        inputsEP.KYMARK = kyhat[ikyhat_mark]

        if (printout)
            println("imark_min <= nfactor ", inputsEP.FACTOR_IN, " ", inputsEP.WIDTH_IN, " ", inputsEP.KYMARK)
        end
    end

    if (printout)
        println("ir: ", inputsEP.IR, " iefwid_mark: ", iefwid_mark, " ikyhat_mark: ", ikyhat_mark)
    end

    # Return to the driver these values and the final growthrate values of the last scan.
    return growthrate, inputsEP, inputsPR

    # This function will be done for however many radii you are testing. These values do not interact in the driver. 
end

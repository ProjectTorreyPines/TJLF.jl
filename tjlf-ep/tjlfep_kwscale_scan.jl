#Temp-calls for working on file:
#=
using MPI
using SparseArrays
include("../tjlf-ep/TJLFEP.jl")
using .TJLFEP
=#
using MPI
"""
kwscale_scan: Each MPI process runs this code for a certain amount of rounds dependent on how
many processes are being ran. 

Unlike the Fortran version where TGLFEP_interface is used, for now I am using the
InputTJLFEP and InputTJLF structs as the combined profile. Otherwise I would need to 
cross-reference between them.

Inputs: inputsEP::InputTJLFEP, inputsPR::InputTJLF

Outputs: N/A yet.

"""
function kwscale_scan(inputsEP::InputTJLFEP{Float64}, inputsPR::profile{Float64}) 
    # The inputsPR will likely be changed to ::profile when it is complete.
    
    # These are for testing purposes:
    #baseDirectory = "/Users/benagnew/TJLF.jl/outputs/tglfep_tests/input.MTGLF"
    #inputsPR = readMTGLF(baseDirectory)
    
    nfactor = 10
    nefwid = 10
    nkyhat = 5
    nkwf = nfactor*nefwid*nkyhat
    k_max = 4

    TJLFEP_COMM = MPI.COMM_WORLD # This will also be modified later.
    #MPI.Init() #This will likely be deleted once this is a function
    np = MPI.Comm_size(TJLFEP_COMM)
    id = MPI.Comm_rank(TJLFEP_COMM)
    # etc. I'll deal this with later

    kyhat_min = 0.0
    kyhat_max = 1.0

    growthrate = 0.0
    growthrate_out = 0.0
    frequency = 0.0
    frequency_out = 0.0

    # Very temporary:
    #baseDirectory = "/Users/benagnew/TJLF.jl/outputs/tglf_regression/tglf01/"
    #inputTJLF = readTGLFEP(baseDirectory*"input.tglf")
    

    # I don't believe I need the _out versions as the the originals don't need
    # to be maintained after all reductions. And Allreduce! will just modify these
    # originals themselves. 
    lkeep_i = fill(true, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    #lkeep_i_out = fill(true, (nkyhat, nefwid, nfactor, inputsPR.NMODES))
    ltearing_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    #ltearing_i_out = fill(false, (nkyhat, nefwid, nfactor, inputsPR.NMODES)) # I don't know if I really need to allocate this in Julia
    l_th_pinch_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_i_pinch_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_e_pinch_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_EP_pinch_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_max_outer_panel_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))
    l_QL_ratio_i = fill(false, (nkyhat, nefwid, nfactor, inputsEP.NMODES))

    f0 = 0.0
    f1 = inputsEP.FACTOR_IN
    w0 = inputsEP.WIDTH_MIN
    w1 = inputsEP.WIDTH_MAX
    kyhat_min = 0.0 # These aren't declared from input
    kyhat_max = 1.0
    kyhat0 = kyhat_min
    kyhat1 = kyhat_max
    factor = fill(NaN, nfactor)
    efwid = fill(NaN, nefwid)
    kyhat = fill(NaN, nkyhat)
    
    #for k = 1:k_max
        # The following three loops establish equidistant spacing in factor, efwid, and kyhat
        # on the range of each we are interested in. 
        for i = 1:nfactor
            factor[i] = (f1-f0)/nfactor*i+f0
        end
        for i = 1:nefwid
            efwid[i] = (w1-w0)/nefwid*i+w0
            println(efwid[i])
        end
        for i = 1:nkyhat
            kyhat[i] = (kyhat1-kyhat0)/nkyhat+kyhat0
        end

        for i = 1+id:np:nkwf-480 # temporarily for speed purposes
            l_wavefunction_out = 0

            # The following 3 statments define each combination of ikyhat, iefwid, and ifactor.
            ikyhat = Int(floor((i-1)/(nefwid*nfactor))+1)
            iefwid = Int(floor(1.0*mod(i-1, nefwid*nfactor)/nfactor)+1)
            ifactor = mod(i-1, nfactor)+1 # this should be 1 for the first ten and then 2 for the next 10 and on.... Okay, this works.

            inputsEP.FACTOR_IN = factor[ifactor]
            inputsEP.KYHAT_IN = kyhat[ikyhat]
            inputsEP.WIDTH_IN = efwid[iefwid]
            #println("WIDTH_IN at pass ", i, ": ", inputsEP.WIDTH_IN)
            #println("np: ", np)

            #println(kyhat_in)
            #println(typeof(kyhat_in))

            #str_sf
            #Then writing the out.wavefunctions

            TJLFEP_ky(inputsEP, inputsPR)

            #Following this is the run of getting:
            #for n = 1:NMODES #Change NMODES later!
                #=
                growthrate(ikyhat, iefwid, ifactor, n) = get_growthrate(n)
                frequency(ikyhat,iefwid,ifactor,n)  = get_frequency(n)
                lkeep_i(ikyhat,iefwid,ifactor,n) = lkeep(n)
                ltearing_i(ikyhat,iefwid,ifactor,n) = ltearing(n)
                l_th_pinch_i(ikyhat,iefwid,ifactor,n) = l_th_pinch(n)
                l_i_pinch_i(ikyhat,iefwid,ifactor,n) = l_i_pinch(n)
                l_e_pinch_i(ikyhat,iefwid,ifactor,n) = l_e_pinch(n)
                l_EP_pinch_i(ikyhat,iefwid,ifactor,n) = l_EP_pinch(n)
                l_max_outer_panel_i(ikyhat,iefwid,ifactor,n) = l_max_outer_panel(n)
                l_QL_ratio_i(ikyhat,iefwid,ifactor,n) = l_QL_ratio(n)
                =#
            #end
            

        end # end of MPI process collection
        
        # Ending here are this point for testing purposes:
        #=

        MPI.Barrier(TGLFEP_COMM)

        MPI.Allreduce!(growthrate, MPI.SUM, TJLFEP_COMM)
        MPI.Allreduce!(frequency, MPI.SUM, TJLFEP_COMM)
        MPI.Allreduce!(lkeep_i, MPI.LAND, TJLFEP_COMM)
        MPI.Allreduce!(ltearing_i, MPI.LOR, TJLFEP_COMM)
        MPI.Allreduce!(l_th_pinch_i, MPI.LOR, TJLFEP_COMM)
        MPI.Allreduce!(l_i_pinch_i, MPI.LOR, TJLFEP_COMM)
        MPI.Allreduce!(l_e_pinch_i, MPI.LOR, TJLFEP_COMM)
        MPI.Allreduce!(l_max_outer_panel_i, MPI.LOR, TJLFEP_COMM)
        MPI.Allreduce!(l_QL_ratio_i, MPI.LOR, TJLFEP_COMM)

        imark = fill(nfactor+1, (nkyhat, nefwid)) 
        for ikyhat = 1:nkyhat
            for iefwid = 1:nefwid
                for ifactor = 1:nfactor
                    for n = 1:NMODES # Change
                        if (lkeep_i[ikyhat, iefwid, ifactor, n])
                            imark[ikyhat, iefwid] = ifactor
                            break
                        end
                    end # n
                    if (imark[ikyhat,iefwid] <= nfactor) break end
                end # ifactor
            end # iefwid
        end #ikyhat

        imark_min = nfactor + 1
        for ikyhat = 1:nkyhat
            for iefwid = 1:nefwid
                imark_min = min(imark[ikyhat, iefwid], imark_min)
            end
        end

        # These parts are relatively easy to translate but less so to
        # both understand easily and run fluidly. There's a lot that
        # needs to be defined by the MPI processes. Running it at one ky
        # would probably help with the single process.
        fmark = 1.0e20
        gmark = 0.0
        f_guess_mark = 1.0e20
        f_guess = fill(NaN, (nkyhat, nefwid))
        gamma_mark_i_1 = fill(NaN, (nkyhat, nefwid))
        gamma_mark_i_2 = fill(NaN, (nkyhat, nefwid))
        f_mark_i = fill(NaN, (nkyhat, nefwid))
        if (imark_min <= nfactor)
            lkeep_ref = fill(false, (nkyhat, nefwid))
            ikyhat_mark = nkyhat
            iefwid_mark = nefwid
            ikyhat_write = floor(nkyhat/2) # not sure the use of these yet
            iefwid_write = floor(nefwid/2)
            for ikyhat = 1:nkyhat
                for iefwid = 1:nefwid
                    imark_ref = nfactor-1
                    if (imark[ikyhat, iefwid] < nfactor)
                        imark_ref = imark[ikyhat, iefwid]
                        lkeep_ref[ikyhat, iefwid] = true
                    else
                        imark_ref = imark[ikyhat, iefwid] - 1
                        if (imark[ikyhat, iefwid] == nfactor)
                            lkeep_ref[ikyhat,iefwid] = true
                        end
                    end
                    f_g1 = factor[imark_ref]
                    f_g2 = factor[imark_ref+1]
                    gamma_g1 = -2.0
                    gamma_g2 = -1.0
                    for n = 1:NMODES # CHANGE LATER!
                        if (lkeep_i[ikyhat, iefwid, imark_ref, n])
                            gamma_g1 = max(maximum(growthrate[ikyhat, iefwid, imark_ref, n]), gamma_g1)
                        end
                        if (lkeep_i[ikyhat, iefwid, imark_ref+1, n])
                            gamma_g2 = max(maximum(growthrate[ikyhat, iefwid, imark_ref+1, n]), gamma_g2)
                        end
                    end
                    f_guess[ikyhat, iefwid] = (gamma_g1*f_g2-gamma_g2*f_g1)/(gamma_g1-gamma_g2)
                    gamma_mark_i_1[ikyhat, iefwid] = gamma_g1
                    gamma_mark_i_2[ikyhat, iefwid] = gamma_g2
                    f_mark_i[ikyhat, iefwid] = f_g1
                end # iefwid
            end # ikyhat
            for ikyhat = 1:nkyhat
                for iefwid = 1:nefwid
                    if ((lkeep_ref[ikyhat, iefwid]) && (gamma_mark_i_1[ikyhat, iefwid] < 0.95*gamma_mark_i_2[ikyhat, iefwid]) &&
                        (f_mark_i[ikyhat, iefwid] < fmark))
                        ikyhat_mark = ikyhat
                        ikyhat_write = ikyhat
                        iefwid_mark = iefwid
                        iefwid_write = iefwid
                        gmark = gamma_mark_i_1[ikyhat, iefwid]
                        fmark = f_mark_i[ikyhat, iefwid]
                        f_guess_mark = f_guess[ikyhat, iefwid]
                    end
                end # iefwid
            end # ikyhat
        end # if ending (mode found)

        # Next is the writing of the out.scalefactor files.
        =#
    #end

        #==#
end

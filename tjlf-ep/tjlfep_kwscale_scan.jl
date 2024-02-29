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
function kwscale_scan(inputsEP::InputTJLFEP{Float64}, inputsPR::profile{Float64}, COMM_IN::MPI.Comm) 
    # The inputsPR will likely be changed to ::profile when it is complete.
    
    # These are for testing purposes:
    #baseDirectory = "/Users/benagnew/TJLF.jl/outputs/tglfep_tests/input.MTGLF"
    #inputsPR = readMTGLF(baseDirectory)
    
    nfactor = 10
    nefwid = 10
    nkyhat = 5
    nkwf = nfactor*nefwid*nkyhat
    k_max = 1
    l_write_out = true

    TJLFEP_COMM = COMM_IN # This will also be modified later.
    #MPI.Init() #This will likely be deleted once this is a function
    np = MPI.Comm_size(TJLFEP_COMM)
    id = MPI.Comm_rank(TJLFEP_COMM)
    # etc. I'll deal this with later

    kyhat_min = 0.0
    kyhat_max = 1.0

    growthrate = zeros(Float64, nkyhat, nefwid, nfactor, inputsEP.NMODES)
    growthrate_out = 0.0
    frequency = zeros(Float64, nkyhat, nefwid, nfactor, inputsEP.NMODES)
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
    
    for k = 1:k_max
        # The following three loops establish equidistant spacing in factor, efwid, and kyhat
        # on the range of each we are interested in. 
        for i = 1:nfactor
            factor[i] = (f1-f0)/nfactor*i+f0
        end
        for i = 1:nefwid
            efwid[i] = (w1-w0)/nefwid*i+w0
            #println(efwid[i])
        end
        for i = 1:nkyhat
            kyhat[i] = (kyhat1-kyhat0)/nkyhat*i+kyhat0
        end

        for i = 1+id:np:nkwf # temporarily for speed purposes # :np:
            l_wavefunction_out = 0

            # The following 3 statments define each combination of ikyhat, iefwid, and ifactor.
            ikyhat = Int(floor((i-1)/(nefwid*nfactor))+1)
            iefwid = Int(floor(1.0*mod(i-1, nefwid*nfactor)/nfactor)+1)
            ifactor = mod(i-1, nfactor)+1 # this should be 1 for the first ten and then 2 for the next 10 and on.... Okay, this works.

            inputsEP.FACTOR_IN = factor[ifactor] # Just used in mapping
            inputsEP.KYHAT_IN = kyhat[ikyhat] # Set equal to ky_in
            inputsEP.WIDTH_IN = efwid[iefwid]
            #kyIndex = ikyhat
            #println("WIDTH_IN at pass ", i, ": ", inputsEP.WIDTH_IN)
            #println("np: ", np)

            #println(kyhat_in)
            #println(typeof(kyhat_in))

            #str_sf
            #Then writing the out.wavefunctions

            gamma_out, freq_out = TJLFEP_ky(inputsEP, inputsPR) # This should be running over just one ky. The problem is that when tjlf_run is called,
            # it is performing the entire operation over all ky. I don't want to run this a bunch of times since it gets rid of the purposes
            # of using MPI. 

            # One of the major things I want to look into tomorrow is the difference between the Fortran and Julia running versions of ky. TJLF and TGLF run differently
            # due to TJLF being able to run all of the ky in parallel while TGLF runs them sequentially. This is where ky_index comes in for TJLF, which is not present in
            # TGLF. I want to figure out how I'm going to need to access the information that I need (see below commented-out) for each kyhat. Presumably in the Fortran version
            # which is directly translated below, the ikyhat is being selected for each i, but the TJLF (check on the veracity of this) runs all ky's in parallel? So do I need
            # to set these in a new fashion rather?

            #Following this is the run of getting:
            for n = 1:inputsEP.NMODES
                growthrate[ikyhat,iefwid,ifactor,n] = gamma_out[n]
                #println(growthrate[ikyhat,iefwid,ifactor,n])
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
            println("Iteration ", i, ", id ", id)
        end # end of MPI process collection
        
        # Ending here are this point for testing purposes:

        MPI.Barrier(TJLFEP_COMM)

        # Allreduce performs the operation (2nd parameter) between each of the processes in a specified color Group
        # determined by TJLFEP_COMM, then distributes in back to each of the processes individually.

        MPI.Allreduce!(growthrate, MPI.SUM, TJLFEP_COMM)
        MPI.Allreduce!(frequency, MPI.SUM, TJLFEP_COMM)
        MPI.Allreduce!(lkeep_i, MPI.LAND, TJLFEP_COMM)
        MPI.Allreduce!(ltearing_i, MPI.LOR, TJLFEP_COMM)
        MPI.Allreduce!(l_th_pinch_i, MPI.LOR, TJLFEP_COMM)
        MPI.Allreduce!(l_i_pinch_i, MPI.LOR, TJLFEP_COMM)
        MPI.Allreduce!(l_e_pinch_i, MPI.LOR, TJLFEP_COMM)
        MPI.Allreduce!(l_max_outer_panel_i, MPI.LOR, TJLFEP_COMM)
        MPI.Allreduce!(l_QL_ratio_i, MPI.LOR, TJLFEP_COMM)

        imark = fill(0, (nkyhat, nefwid))
        imark .= nfactor + 1
        for ikyhat = 1:nkyhat
            for iefwid = 1:nefwid
                for ifactor = 1:nfactor
                    for n = 1:inputsEP.NMODES # For a specified mode, 
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
        lkeep_ref = fill(false, (nkyhat, nefwid))
        if (imark_min <= nfactor)
            
            ikyhat_mark = nkyhat
            iefwid_mark = nefwid
            ikyhat_write = floor(Int, nkyhat/2) # not sure the use of these yet
            iefwid_write = floor(Int, nefwid/2)
            for ikyhat = 1:nkyhat
                for iefwid = 1:nefwid
                    lkeep_ref[ikyhat, iefwid] = false
                    imark_ref = nfactor-1
                    #println("pass 1: ", imark_ref)
                    if (imark[ikyhat, iefwid] < nfactor)
                        imark_ref = imark[ikyhat, iefwid]
                        #println("route 1: ", imark_ref)
                        lkeep_ref[ikyhat, iefwid] = true
                    else
                        imark_ref = imark[ikyhat, iefwid] - 1
                        #println("route 2: ", imark_ref)
                        if (imark[ikyhat, iefwid] == nfactor)
                            lkeep_ref[ikyhat,iefwid] = true
                        end
                    end
                    #println("end: ", imark_ref)
                    f_g1 = factor[imark_ref]
                    if (imark_ref+1 < 11)
                        f_g2 = factor[imark_ref+1]
                    else
                        f_g2 = 0
                    end
                    gamma_g1 = -2.0
                    gamma_g2 = -1.0
                    for n = 1:inputsEP.NMODES
                        if (lkeep_i[ikyhat, iefwid, imark_ref, n])
                            gamma_g1 = max(maximum(growthrate[ikyhat, iefwid, imark_ref, n]), gamma_g1)
                        end
                        if (imark_ref+1 < 11 && lkeep_i[ikyhat, iefwid, imark_ref+1, n])
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
        
        # These will be produced whenever I run from a specified directory, say
        # via a batch file or some execution command for TJLF-EP. Running by REPL
        # does create these files, but it places them in the TJLF.JL directory,
        # which isn't where I will want to store them in the end.
        if (l_write_out && id == 0)
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
        end
    end

    return growthrate
end

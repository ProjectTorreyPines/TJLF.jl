# More temp calls:
#include("../tjlf-ep/TJLFEP.jl")
#using .TJLFEP
#using .TJLFEP: convert_input
#using MPI

function TJLFEP_ky(inputsEP::Options{Float64}, inputsPR::profile{Float64}, str_wf_file::String, l_wavefunction_out::Int, printout::Bool = true) #, factor_in::Int64, kyhat_in::Int64, width_in::Int64)

    # Temp Defs:
    
    #color = 0
    #kyhat_in = 3
    # Temp Struct Inputs:
    #filename = "/Users/benagnew/TJLF.jl/outputs/tglfep_tests/input.MTGLF"
    #temp = readMTGLF(filename)
    #inputsPR = temp[1]
    #irexp2 = temp[2]
    #filename = "/Users/benagnew/TJLF.jl/outputs/tglfep_tests/input.TGLFEP"
    #inputsEP = readTGLFEP(filename, irexp2)
    #inputsEP.IR = inputsEP.IR_EXP[color+1]
    #inputsEP.MODE_IN = 2
    #inputsEP.KY_MODEL = 3

    #========================================#

    inputTJLF = TJLF_map(inputsEP, inputsPR)

    #println("GradBFactor:")
    #println(inputTJLF.GRADB_FACTOR)

    inputTJLF.USE_TRANSPORT_MODEL = false # this will just run a single ky
    #println(inputTJLF.USE_TRANSPORT_MODEL)
    if (inputsEP.PROCESS_IN == 0)
        # See TGLF-EP if needed.
    end
    
    inputTJLF.KYGRID_MODEL = 0

    inputTJLF.NMODES = inputsEP.NMODES

    inputTJLF.NBASIS_MIN = inputsEP.N_BASIS
    inputTJLF.NBASIS_MAX = inputsEP.N_BASIS

    inputTJLF.NXGRID = 32

    inputTJLF.WIDTH = inputsEP.WIDTH_IN
    inputTJLF.FIND_WIDTH = false



    # Corrections for TJLF specifically: (see main.jl after mainsub call)

    inputTJLF.USE_AVE_ION_GRID = false
    inputTJLF.WIDTH_SPECTRUM .= inputTJLF.WIDTH # see tjlf_read_input.jl from TJLF.jl.
    #println("WIDTH_SPECTRUM: ", inputTJLF.WIDTH_SPECTRUM)
    #println("WIDTH: ", inputTJLF.WIDTH)
    inputTJLF.FIND_EIGEN = true # in all inputs for tjlf this is set to true
    inputTJLF.RLNP_CUTOFF = 18.0 # in all inputs for tjlf this is set to 18.0
    inputTJLF.BETA_LOC = 0.0 # This one I am very unsure of. Some 0.0, some 1.0. 
    inputTJLF.DAMP_PSI = 0.0 # in all inputs for tjlf this is set to 0.0
    inputTJLF.DAMP_SIG = 0.0 # in all inputs for tjlf this is set to 0.0
    inputTJLF.WDIA_TRAPPED = 0.0

   
    # n_out::Int
    # EP_QL_e_flux::Float32
    # ef_phi_norm::Float32

    if inputTJLF.SAT_RULE == 2 || inputTJLF.SAT_RULE == 3 # From read_input, which is skipped over in this path of running TJLF
        inputTJLF.UNITS = "CGYRO"
        ### WTF
        inputTJLF.XNU_MODEL = 3
        inputTJLF.WDIA_TRAPPED = 1.0
    end

    inputTJLF.KX0_LOC = 0.0
    #println("inputTJLF: ")
    #println(typeof(inputTJLF))

    convInput = convert_input(inputTJLF, inputTJLF.NS, inputTJLF.NKY)
    #println("convInput: ")
    #println(typeof(convInput))


    # A few test print statements which most likely aren't needed anymore:
    #=
    if (inputsEP.IR == 101)
        io78 = open("run.out.second", "a")
        fieldnames_in = fieldnames(typeof(convInput))
        for field in fieldnames_in
            println(io78, getfield(convInput, field), ": ", field)
        end
        close(io78)
    end
    =#

    #=if (inputsEP.IR == 2)
        println(inputTJLF.KX0_LOC)
        println(convInput.KX0_LOC)
        println("====")
    end


    if (inputsEP.IR == 2)
        io77 = open("run.out.second", "a")
        println(io77, convInput)
        println(io77, convInput.BETA_LOC)
        println(io77, convInput.KX0_LOC)
        println(io77, convInput.PSI)
        close(io77)
    end

    if (inputsEP.IR == 2)
        println("Second")
        println(inputTJLF.KX0_LOC)
        println(convInput.KX0_LOC)
        println("====")
    end=#

    #if (inputsEP.IR == 2)
    #    println(convInput.IBRANCH)
    #    println(convInput.WIDTH)
    #end

    # Run TJLF and return QLweight and eigenvalues:
    gamma_out, freq_out, particle_QL_out, energy_QL_out, stress_par_QL_out, exchange_QL_out, field_weight_out, satParams, nmodes_out = TJLF.run(convInput)




    #println(typeof(field_weight_out))
    #println(typeof(satParams.y))


    inputTJLF = revert_input(convInput, convInput.NS, convInput.NKY)
    


    # Next is the get_growthrate stuff. I now need to make sure that run_TJLF is giving me all the information I need to continue this.

    g = fill(NaN, inputTJLF.NMODES)
    f = fill(NaN, inputTJLF.NMODES)
    for n = 1:inputTJLF.NMODES
        g[n] = gamma_out[n]
        f[n] = freq_out[n]
    end

    # println("before ", gamma_out)
    #GAMMA STILL NON-ZERO UP TO HERE

    # Establishes the lkeep vector. LKEEP is defaulted to all true as nothing is rejected yet.
    # This states that if the frequency that came from the converted test of TJLF is less than
    # the cutoff, the mode is kept. It also then requires that the growthrate is larger than 1e-7.
    inputsEP.LKEEP .= true

    for n = 1:inputTJLF.NMODES
        inputsEP.LKEEP[n] = (f[n] < inputsEP.FREQ_AE_UPPER)
        inputsEP.LKEEP[n] = (inputsEP.LKEEP[n] && (g[n] > inputsEP.GAMMA_THRESH))
    end


    if (inputsEP.IR == 101 && printout)

        println("----------")
        
        #println(satParams.y)
        #println(satParams.theta)
        #println("----------")
        #println("nmodes_out: ", nmodes_out)
        #println("nbasis: ",inputsEP.N_BASIS)
        #println("==========")

    end

    # This function was translated within TJLF so as to get the wavefunction.
    ms = 128
    max_plot = Int(18*ms/8+1)
    wavefunction, angle, nplot, nmodes, nmodes_out = TJLF.get_wavefunction(convInput, satParams, field_weight_out, nmodes_out)


    inputsEP.LTEARING .= false
    inputsEP.L_I_PINCH .= false
    inputsEP.L_E_PINCH .= false
    inputsEP.L_TH_PINCH .= false
    inputsEP.L_QL_RATIO .= false
    # inputsEP.L_MAX_OUTER_PANEL .= false
    x_tear_test::Vector{Float64} = fill(0.0, 4)
    abswavefunction = abs.(wavefunction)

    

    #absdiffwavefunction = similar(wavefunction)
    ms = 128
    npi = 9
    np = Int(ms/8)
    nb = inputTJLF.NBASIS_MAX
    igeo = 1 # Hard-Coded for now as with the previous LS functions:
    max_plot = Int(18*ms/8+1) # 289 length vector
    maxmodes = inputTJLF.NMODES

    i_QL_cond_flux = fill(0.0, 4)
    e_QL_cond_flux = fill(0.0, 4)
    QL_flux_ratio = fill(0.0, 4)
    EP_conv_frac = fill(0.0, 4)
    theta_2_moment = fill(0.0, 4)
    # TGLFEP hard-codes nmodes = 4, so that is why these are all defined like this.
    DEP = fill(NaN, 4)
    chi_th = fill(NaN, 4)
    chi_i = fill(NaN, 4)
    chi_i_cond = fill(NaN, 4)
    chi_e = fill(NaN, 4)
    chi_e_cond = fill(NaN, 4)

    #println(wavefunction)
    #=if (inputsEP.IR == 101)
        println("===nmodes_out===")
        println(nmodes_out)
    end=#

    for n = 1:inputsEP.NMODES
        # The use of NMODES here is slightly confusing but needed as the Fortran uses assumed values for modes
        # which did not satisfy the criteria (nmodes_out). This means that the loop must continue for
        # those past nmodes_out so as to be consistent with the flags. This is why n <= nmodes_out is used
        # multiple times in this loop.

        #nul = abswavefunction
        wave_max = maximum(abs.(wavefunction[n,1,:]))+1.0E-3
        
        wave_max_loc = argmax(abs.(wavefunction[n,1,:]))
        n_balloon_pi = floor(Int, (max_plot-1)/9) # 32
        i_mid_plot = floor(Int, (max_plot-1)/2+1) # 145
        # inputsEP.L_MAX_OUTER_PANEL[n] = (wave_max_loc < (i_mid_plot-n_balloon_pi)) || (wave_max_loc > (i_mid_plot+n_balloon_pi))
        # So long as wave_max_loc is between 113 and 177, this isn't rejected for this.
        theta_2_moment[n] =0.0
        ef_phi_norm = 0.0

        if (inputsEP.IR == 101 && n == 4 && printout)
            println("wave_max: ", wave_max)
            println("wave_max_loc: ", wave_max_loc)
            println("mode: ", n)
            #println("inputsEP.L_MAX_OUTER_PANEL[n]: ", inputsEP.L_MAX_OUTER_PANEL[n])
            println("-----")
            println("wavefunction: ", abs.(wavefunction[n,1,wave_max_loc:wave_max_loc+10]))
            println("=====")
        end


        if (n <= nmodes_out)#used to be nmodes_out
            for i = 1:max_plot # Finding the maximum value of this abs value of difference div wave_max
                absdiffwavefunction::Float64 = abs.(wavefunction[n,1,i]-wavefunction[n,1,max_plot+1-i])
                x_tear_test[n] = max(x_tear_test[n], absdiffwavefunction/wave_max)
                ef_phi_norm += abs(wavefunction[n,1,i])
                theta_2_moment[n] += (9 * π * (-1.0 + (2.0 * (i - 1)) / (max_plot - 1)))^2 * abs(wavefunction[n, 1, i])
                # println("xteartest= ", x_tear_test)
                # println("theta_2_moment = ", theta_2_moment) 
                # println("ef_phi_norm = ", ef_phi_norm)  
            end
            theta_2_moment[n] /= ef_phi_norm

            if (x_tear_test[n] > 1.0E-1)
                inputsEP.LTEARING[n] = true
                
            end
        end 
        EP_QL_e_flux = 0.0
        EP_QL_flux = 0.0
        i_QL_flux = 0.0
        i_QL_cond_flux[n] = 0.0
        i_eff_grad = 0.0
        e_QL_flux = 0.0
        e_QL_cond_flux[n] = 0.0
        th_QL_flux = 0.0
        th_eff_grad = 0.0

        if (n <= nmodes_out)
            for jfields = 1:3
                # println("size(particle_QL_out ", size(particle_QL_out))
                EP_QL_flux = EP_QL_flux + particle_QL_out[jfields, inputsEP.IS_EP + 1, n]
                # io97 = open("radius.test.out", "a")
                # println(io97, "n = ", n)
                # println(io97, "jfields = ", jfields)
                # println(io97, "particle_QL_out[jfields, inputsEP.IS_EP + 1, n]= ", particle_QL_out[jfields, inputsEP.IS_EP + 1, n])
                # println(io97, "energy_QL_out[jfields, inputsEP.IS_EP + 1, n]= ", energy_QL_out[jfields, inputsEP.IS_EP + 1, n])
                # close(io97)
                EP_QL_e_flux = EP_QL_e_flux + energy_QL_out[jfields, inputsEP.IS_EP + 1, n]
                e_QL_flux = e_QL_flux + energy_QL_out[jfields, 1, n]
                e_QL_cond_flux[n] = e_QL_cond_flux[n] + energy_QL_out[jfields, 1, n] - 1.5*inputTJLF.TAUS[1]*particle_QL_out[jfields, 1, n]
                for j_ion = 2:inputsEP.IS_EP
                    i_QL_flux = i_QL_flux + energy_QL_out[jfields, j_ion, n]
                    i_QL_cond_flux[n] = i_QL_cond_flux[n] + energy_QL_out[jfields, j_ion, n] - 1.5*inputTJLF.TAUS[j_ion]*particle_QL_out[jfields, j_ion, n]
                end         
            end
        end

        for j_ion = 2:inputsEP.IS_EP
            i_eff_grad = i_eff_grad + inputTJLF.TAUS[j_ion]*inputTJLF.AS[j_ion]*inputTJLF.RLTS[j_ion]
        end

        th_eff_grad = i_eff_grad + inputTJLF.RLTS[1]*inputTJLF.AS[1]
        th_QL_flux = i_QL_cond_flux[n] + e_QL_cond_flux[n]
        DEP[n] = EP_QL_flux / (inputTJLF.AS[inputsEP.IS_EP+1]*inputTJLF.RLNS[inputsEP.IS_EP+1])
        chi_e[n] = e_QL_flux / (inputTJLF.RLTS[1]*inputTJLF.AS[1])
        chi_e_cond[n] = e_QL_cond_flux[n] / (inputTJLF.RLTS[1]*inputTJLF.AS[1])
        chi_i[n] = i_QL_flux / i_eff_grad
        chi_i_cond[n] = i_QL_cond_flux[n] / i_eff_grad
        chi_th[n] = th_QL_flux / th_eff_grad

        if (n <= nmodes_out) 
                # QL_flux_ratio[n] = (EP_QL_flux/inputTJLF.AS[inputsEP.IS_EP+1])/(abs(i_QL_cond_flux[n])/(inputTJLF.AS[1]-inputTJLF.AS[inputsEP.IS_EP+1]))
            QL_flux_ratio[n] = EP_QL_e_flux / abs(i_QL_cond_flux[n])
            EP_conv_frac[n] = EP_QL_flux * 1.5 * inputTJLF.TAUS[inputsEP.IS_EP+1] / EP_QL_e_flux
        else
            QL_flux_ratio[n] = 1e20
        end

        if (chi_i[n] < 0.0) inputsEP.L_I_PINCH[n] = true end
        if (chi_e[n] < 0.0) inputsEP.L_E_PINCH[n] = true end
        if (chi_th[n] < 0.0) inputsEP.L_TH_PINCH[n] = true end
        if (DEP[n] < 0.0) inputsEP.L_EP_PINCH[n] = true end
        if (QL_flux_ratio[n] < inputsEP.QL_RATIO_THRESH) inputsEP.L_QL_RATIO[n] = true end
        if (theta_2_moment[n] > inputsEP.THETA_SQ_THRESH) inputsEP.L_THETA_SQ[n] = true end
        
       
    end

    #=if (testid)
        println(chi_th)
    end=#
    
    # Runs over all modes (4) for the lkeep vector and checks if each flag is false. If lkeep false, all will be false. 
    # If lkeep is true, if the flag is false, lkeep will stay true; if lkeep is true and the flag is true, lkeep will be turned back to false
    # This essentiall means that if the rejection flag is turned on, anything marked for rejection will be rejected.
    for n = 1:inputTJLF.NMODES
        if (inputsEP.REJECT_TEARING_FLAG == 1) inputsEP.LKEEP[n] = (inputsEP.LKEEP[n] && !inputsEP.LTEARING[n]) end
        if (inputsEP.REJECT_I_PINCH_FLAG == 1) inputsEP.LKEEP[n] = (inputsEP.LKEEP[n] && !inputsEP.L_I_PINCH[n]) end
        if (inputsEP.REJECT_E_PINCH_FLAG == 1) inputsEP.LKEEP[n] = (inputsEP.LKEEP[n] && !inputsEP.L_E_PINCH[n]) end
        if (inputsEP.REJECT_TH_PINCH_FLAG == 1) inputsEP.LKEEP[n] = (inputsEP.LKEEP[n] && !inputsEP.L_TH_PINCH[n]) end
        if (inputsEP.REJECT_EP_PINCH_FLAG == 1) inputsEP.LKEEP[n] = (inputsEP.LKEEP[n] && !inputsEP.L_EP_PINCH[n]) end
        # if (inputsEP.ROTATIONAL_SUPPRESSION_FLAG == 1) inputsEP.LKEEP[n] = (inputsEP.LKEEP[n] && !inputsEP.L_MAX_OUTER_PANEL[n]) end
        inputsEP.LKEEP[n] = (inputsEP.LKEEP[n] && !inputsEP.L_QL_RATIO[n])
        inputsEP.LKEEP[n] = (inputsEP.LKEEP[n] && !inputsEP.L_THETA_SQ[n])
    end
    #extra print to test variables

    # println("absdiffwavefunction",absdiffwavefunction)
    


    # Next is writing the wavefunction files themselves:
    if (l_wavefunction_out == 1) # nplot = max_plot_out; nfields = 1 by def.
        io6 = open(str_wf_file, "w")
        println(io6, "nmodes=", nmodes_out, "TGLF INPUT NMODES = ",inputsEP.NMODES, " nfields=", " max_plot=", nplot)
        println(io6, "ky=", inputTJLF.KY, " width=", inputTJLF.WIDTH)
        println(io6, "theta     ", "((Re(field_i), Im(field_i),i=(1,nfields)),j=1,nmodes)")
        println(io6, "Tearing metric: ", x_tear_test)
        println(io6, "DEP: ", DEP)
        println(io6, "chi_th: ", chi_th)
        println(io6, "chi_i: ", chi_i)
        println(io6, "chi_i_cond: ", chi_i_cond)
        println(io6, "chi_e: ", chi_e)
        println(io6, "chi_e_cond: ", chi_e_cond)
        println(io6, "i_QL_cond_flux: ", i_QL_cond_flux)
        println(io6, "e_QL_cond_flux: ", e_QL_cond_flux)
        println(io6, "QL_ratio: ", QL_flux_ratio)
        println(io6, "EP QL convection fracton: ", EP_conv_frac)
        println(io6, "<theta^2>: ", theta_2_moment)
        println(io6, "lkeep: ", inputsEP.LKEEP)
        # Renormalize and adjust phases:
        n_out = 0
        for n = 1:inputsEP.NMODES
            max_phi = maximum(abswavefunction[n,1,:])
            max_apar = maximum(abswavefunction[n,2,:])
            max_field = maximum([max_phi, max_apar])
            phase = atan(imag(wavefunction[n,1,Int((nplot+1)/2)]),real(wavefunction[n,1,Int((nplot+1)/2)]))
            for jfields = 1:2
                z = 0+1im
                wavefunction[n,jfields,:] .= wavefunction[n,jfields,:]/(max_field*exp(z*phase))
            end
            if n_out == 0 && inputsEP.LKEEP[n]
                n_out = n
            end
        end
        #println("LKEEP ", inputsEP.LKEEP)
        #println("n_out = " ,n_out)
        if n_out == 0
            n_out = 1
            println(io6, "No kept modes at nominal write parameters. Showing leading mode.")
        end
        #Write renormalized, re-phased eigenfunctions out to str_wf_file.
        for i = 1:max_plot
            #This doesn't do anything in TGLFEP (???)
            #kcomp = 0
            #wave_write = fill(0.0, 20)
            #for n = 1:nmodes_out
            #    for jfields = 1:3
            #        kcomp += 1
            #        wave_write[kcomp] = real(wavefunction[n,jfields,i])
            #        kcomp += 1
            #        wave_write[kcomp] = imag(wavefunction[n,jfields,i])
            #    end
            #end
            println(io6, angle[i], " ", real(wavefunction[n_out,1,i]), " ", imag(wavefunction[n_out,1,i]), " ",
                    real(wavefunction[n_out,2,i]), " ", imag(wavefunction[n_out,2,i]))
        end
        close(io6)
    end      
    # This is the end of ky.jl. Returning values will likely need to be changed later.
    return gamma_out, freq_out, inputTJLF
end

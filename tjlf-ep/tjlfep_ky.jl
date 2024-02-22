# More temp calls yay:
#include("../tjlf-ep/TJLFEP.jl")
#using .TJLFEP
#using .TJLFEP: convert_input

function TJLFEP_ky(inputsEP::InputTJLFEP{Float64}, inputsPR::profile{Float64}) #, factor_in::Int64, kyhat_in::Int64, width_in::Int64)

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

    if (inputsEP.PROCESS_IN == 0)
        # I probably won't use this process_in value
    end
    
    inputTJLF.KYGRID_MODEL = 0
    # the ky_in is set more fluidly here. Setting this part is mostly for if you're using process_in == 0

    inputTJLF.NMODES = inputsEP.NMODES

    inputTJLF.NBASIS_MIN = inputsEP.N_BASIS
    inputTJLF.NBASIS_MAX = inputsEP.N_BASIS

    inputTJLF.NXGRID = 32

    inputTJLF.WIDTH = inputsEP.WIDTH_IN
    inputTJLF.FIND_WIDTH = false
    #inputTJLF.USE_AVE_ION_GRID = 
    #a = inputTJLF
    # Next is the TJLF call
    # I legit have no clue why this is so finnicky, making sure all the represented functions are available :(
    #typeof(inputTJLF)

    #include("../src/TJLF.jl")

    # Corrections for TJLF specifically: (see main.jl after mainsub call)
    inputTJLF.USE_AVE_ION_GRID = false
    inputTJLF.WIDTH_SPECTRUM .= inputTJLF.WIDTH # see tjlf_read_input.jl from TJLF.jl.
    #println("WIDTH_SPECTRUM: ", inputTJLF.WIDTH_SPECTRUM)
    #println("WIDTH: ", inputTJLF.WIDTH)
    inputTJLF.FIND_EIGEN = true # in all inputs for tjlf this is set to true
    inputTJLF.RLNP_CUTOFF = 18.0 # in all inputs for tjlf this is set to 18.0
    inputTJLF.BETA_LOC = 1.0 # This one I am very unsure of. Some 0.0, some 1.0. 
    inputTJLF.DAMP_PSI = 0.0 # in all inputs for tjlf this is set to 0.0
    inputTJLF.DAMP_SIG = 0.0 # in all inputs for tjlf this is set to 0.0
    inputTJLF.WDIA_TRAPPED = 0.0
    if inputTJLF.SAT_RULE == 2 || inputTJLF.SAT_RULE == 3 # From read_input, which is skipped over in this path of running TJLF
        inputTJLF.UNITS = "CGYRO"
        ### WTF
        inputTJLF.XNU_MODEL = 3
        inputTJLF.WDIA_TRAPPED = 1.0
    end

    #println("inputTJLF: ")
    #println(typeof(inputTJLF))
    convInput = convert_input(inputTJLF, inputTJLF.NS, inputTJLF.NKY)
    #println("convInput: ")
    #println(typeof(convInput))
    #run_TJLF requires a using .TJLF statment...
    #This does in fact make sense. The problem is, I don't know if this will work until
    #I test it from main.jl because that is where I will be running everything. Not from tjlfep_ky.jl or the TWO in-between files kwscale_scan and mainsub.
    #mainsub isn't really a problem but kwscale_scan is because I haven't done much in there yet. I can try to get it to call everything up to the ky scan I suppose.
    QL_weights, eigenvalue, QL_flux_out, flux_spectrum, field_weight_out, satParams = run_TJLF(convInput)
    println(typeof(field_weight_out))
    println(typeof(satParams.y))

    # Next is the get_growthrate stuff. I now need to make sure that run_TJLF is giving me all the information I need to continue this.

    # OK, yes. The output for run_TJLF contains already the frequency and growthrate info I want.
    g = fill(NaN, inputTJLF.NMODES)
    f = fill(NaN, inputTJLF.NMODES)
    for n = 1:inputTJLF.NMODES
        g[n] = eigenvalue[1,n,1]
        f[n] = eigenvalue[1,n,2]
    end

    for n = 1:inputTJLF.NMODES
        inputsEP.LKEEP[n] = f[n] < inputsEP.FREQ_AE_UPPER
        inputsEP.LKEEP[n] = inputsEP.LKEEP[n] && (g[n] > 1.0e-7)
    end

    # I tried my best... Let's see if it was good enough for a first try lol. Probably not but I'mma see.
    wavefunction, angle, nplot, nmodes = get_wavefunction(convInput, satParams, field_weight_out)

    inputsEP.LTEARING .= false
    inputsEP.L_I_PINCH .= false
    inputsEP.L_E_PINCH .= false
    inputsEP.L_TH_PINCH .= false
    inputsEP.L_QL_RATIO .= false
    inputsEP.L_MAX_OUTER_PANEL .= false
    x_tear_test = fill(0, 4)
    abswavefunction = similar(wavefunction)
    abswavefunction .= abs(wavefunction)
    absdifffunction = similar(wavefunction)
    ms = 128
    npi = 9
    np = Int(ms/8)
    nb = inputs.NBASIS_MAX
    igeo = 1 # Hard-Coded for now as with the previous LS functions:
    max_plot = Int(18*ms/8+1) # 289 length vector
    maxmodes = inputs.NMODES

    i_QL_cond_flux = fill(NaN, 4)
    e_QL_cond_flux = fill(NaN, 4)
    QL_flux_ratio = fill(NaN, 4)

    DEP = fill(NaN, 4)
    chi_th = fill(NaN, 4)
    chi_i = fill(NaN, 4)
    chi_i_cond = fill(NaN, 4)
    chi_e = fill(NaN, 4)
    chi_e_cond = fill(NaN, 4)

    for n = 1:inputTJLF.NMODES
        wave_max = maximum(abswavefunction[n,1,:])+1.0e-3
        wave_max_loc = argmax(abswavefunction[n,1,:])
        n_balloon_pi = Int((max_plot-1)/9)
        i_mid_plot = Int(max_plot-1)/2+1
        inputsEP.L_MAX_OUTER_PANEL[n] = (wave_max_loc[1] < (i_mid_plot-n_balloon_pi)) || (wave_max_loc[1] > (i_mid_plo+n_balloon_pi))
        for i = 1:max_plot # Finding the maximum value of this abs value of difference div wave_max
            x_tear_test[n] = max(x_tear_test[n], abs(wavefunction[n,1,i]-wavefunction[n,1,max_plot+1-i])/wave_max)
        end
        if (x_tear_test[n] > 1.0e-1)
            inputsEP.LTEARING[n] = true
        end
        EP_QL_flux = 0.0
        i_QL_flux = 0.0
        i_QL_cond_flux[n] = 0.0
        i_eff_grad = 0.0
        e_QL_flux = 0.0
        e_QL_cond_flux[n] = 0.0
        th_QL_flux = 0.0
        th_eff_grad = 0.0
        for jfields = 1:3 # There's a BUNCH of other functions here which are needed. I don't think they'll be as tedious however. 
            EP_QL_flux = EP_QL_flux #+ get_QL_particle_flux(n, is_EP+1, jfields)
            e_QL_flux = e_QL_flux #+ get_QL_energy_flux(n, 1, jfields)
            e_QL_cond_flux[n] = e_QL_cond_flux[n] #+ get_QL_energy_flux(n, 1, jfields) - 1.5*inputTJLF.TAUS[1]*get_QL_particle_flux(n, 1, jfields)
            for j_ion = 2:inputsEP.IS_EP
                i_QL_flux = i_QL_flux #+ get_QL_energy_flux(n, j_ion, jfields)
                i_QL_cond_flux[n] = i_QL_cond_flux[n] #+ get_QL_energy_flux(n, j_ion, jfields) - 1.5*inputTJLF.TAUS[j_ion]*get_QL_particle_flux(n, j_ion, jfields)
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
        QL_flux_ratio[n] = (EP_QL_flux/inputTJLF.AS[inputsEP.IS_EP+1])/(abs(i_QL_cond_flux[n])/(inputTJLF.AS[1]-inputTJLF.AS[inputsEP.IS_EP+1]))

        if (chi_i[n] < 0.0) inputsEP.L_I_PINCH[n] = true end
        if (chi_e[n] < 0.0) inputsEP.L_E_PINCH[n] = true end
        if (chi_th[n] < 0.0) inputsEP.L_TH_PINCH[n] = true end
        if (DEP[n] < 0.0) inputsEP.L_EP_PINCH[n] = true end
        if (QL_flux_ratio[n] < inputsEP.QL_THRESH_RATIO) inputsEP.L_QL_RATIO[n] = true end
    end

    for n = 1:inputTJLF.NMODES
        if (inputsEP.REJECT_TEARING_FLAG == 1) end
        if (inputsEP.REJECT_I_PINCH_FLAG == 1) end
        if (inputsEP.REJECT_E_PINCH_FLAG == 1) end
        if (inputsEP.REJECT_TH_PINCH_FLAG == 1) end
        if (inputsEP.REJECT_EP_PINCH_FLAG == 1) end
        if (inputsEP.REJECT_MAX_OUTER_PANEL_FLAG == 1) end
        
    end






end
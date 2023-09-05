#
#     parameters: 
#     ky_mix - array of ky (mode number)
#     gamma_mix - array of gamma (net growth rate)
#     kwargs... - dictionary of global variables (rho_ion, alpha_zf_in, grad_r0_out)
#     
#     outputs:
#     vzf_mix - zonal flow mixing rate
#     kymax_mix - ky value at the calculated vzf_mix
#     jmax_mix - index of ky_mix array where vzf_mix is calculated
#
#     finds the maximum of gamma/ky spectrum at low-k values by going through the ky_mix
#     array, then after finding the max, interpolate the value to improve accuracy
#   

##### define the parameter data types?
function get_zonal_mixing(ky_mix,gamma_mix;kwargs...)


    # uses variables from tjlf_dimensions, tglf_global, tglf_species
    # create dictionary
    kwargs=Dict(kwargs)


#
# find the local maximum of gamma_mix/ky_mix with the largest gamma_mix/ky_mix^2
#
    ### kycut is defined with global variables in the python version, but rho_ion is a globalVar
    kycut=0.8/kwargs["rho_ion"]
    kymin = 0.0  

    testmax = 0.0
    jmax_mix = 1
      

    
    kwargs["alpha_zf_in"] < 0.0 ? kymin = 0.173 * √(2.0) / kwargs["rho_ion"] :
    # saturation rules
    if sat_rule_in==2 || sat_rule_in==3
        kycut = kwargs["grad_r0_out"] * kycut
        kymin = kwargs["grad_r0_out"] * kymin
    end

    # find the low and high ky peaks of gamma/ky

    # go through ky_mix except for last value
    # if next val is >= kymin and curr val is less than kycut save it
    # update testmax and max index if necessary

    ##### What should i initialize j1 to???
    j1 = 0
    for j in range(1, length(ky_mix)-1)
        if((ky_mix[j+1] >= kymin) && (ky_mix[j] <= kycut))
            # save index in case no max
            j1=j
            testmax1 = gamma_mix[j]/ky_mix[j]
            # find maxes
            if(testmax1 > testmax)
                testmax = testmax1
                jmax_mix = j
            end
        end
    end

    ###### What do you do if jmax_mix = NOTHING (0 in this case)
    # if testmax is not updated a single time
    testmax==0.0 ? jmax_mix=j1 :

    # no unstable modes in range set kymax index to end of range
    # this is cut of at j1 since a maximum may not exist in the low-k range
    kymax1 = ky_mix[jmax_mix]
    gammamax1 = gamma_mix[jmax_mix]




    # linearly interpolate to find a more accurate low-k maximum gamma/ky
    if(kymax1<kymin)
        kymax1 = kymin
        gammamax1 = gamma_mix[1]
                                +(gamma_mix[2]-gamma_mix[1])*(kymin-ky_mix[1])/(ky_mix[2]-ky_mix[1])
    end

    # determine kymax1 and gammamax1 bounded by the tree points f0,f1,f2
    # use a quadratic fit: f = a + b x + c x^2    to f = gamma/ky centered at jmax1
    # scale it to be quadratic where x goes from 0 to 1
    ##### Are we not worried that jmax_mix is at the edge?
    if(jmax_mix>1 && jmax_mix < j1)
        jmax1 = jmax_mix
        f0 = gamma_mix[jmax1-1]/ky_mix[jmax1-1]
        f1 = gamma_mix[jmax1]/ky_mix[jmax1]
        f2 = gamma_mix[jmax1+1]/ky_mix[jmax1+1]
        deltaky = ky_mix[jmax1+1]-ky_mix[jmax1-1]
        x1 = (ky_mix[jmax1]-ky_mix[jmax1-1])/deltaky
        a = f0
        b = (f1 - f0*(1-x1*x1)-f2*x1*x1)/(x1-x1*x1)
        c = f2 - f0 - b

        # if f0>f1 then f1 is not a local maximum
        if(f0 > f1)
            kymax1 = ky_mix[jmax1-1]
            gammamax1 = f0*kymax1

            #interpolate to find the value of gammamax1 at kymin
            if(kymax1 < kymin)
                kymax1 = kymin
                xmin = (kymin - ky_mix[jmax1-1])/deltaky
                gammamax1 = (a + b*xmin + c*xmin*xmin)*kymin
            end
        end

        # if f0<f1 then f1>f2 due to the maximum search
        # use the quadratic fit to refine the local maximum:
        if(f0<f1)
            xmax = -b/(2.0*c)
            xmin = 0.0

            if(ky_mix[jmax1-1]<kymin)
                xmin = (kymin - ky_mix[jmax1-1])/deltaky
            end

            # if xmax >= 1    use f2 as the maximum
            if(xmax >= 1.0)
                kymax1 = ky_mix[jmax1+1]
                gammamax1 = f2*kymax1
            elseif(xmax <= xmin)   
                # use the quadratic fit to determine gammamax1 at kymin 
                if(xmin > 0.0)
                    kymax1 = kymin
                    gammamax1 = (a + b*xmin + c*xmin*xmin)*kymin
                # if xmax<=0 use f0 as the maximum
                elseif(xmax < 0.0)
                    kymax1 = ky_mix[jmax1-1]
                    gammamax1 = f0*kymax1
                end

            # the conditions f0<f1<f2 and xmin<xmax<1 are satisfied
            # use the quadratic fit to determine gammamax1 and kymax1
            else
                kymax1 = ky_mix[jmax1-1]+deltaky*xmax
                gammamax1 = (a+b*xmax+c*xmax*xmax)*kymax1
            end #xmax tests
        
        end #f0 < f1
    
    end    # jmax_mix > 1

    vzf_mix = gammamax1/kymax1
    kymax_mix = kymax1
    ### commented out in original Fortran code
#            jmax_mix = jmax1

    return vzf_mix, kymax_mix, jmax_mix
 
end
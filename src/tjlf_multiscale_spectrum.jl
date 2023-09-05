

include("tjlf_modules.jl")

# nmix was in parameters, but not used at all?
# vzf_mix,kymax_mix,jmax_mix are output parameters?
function get_zonal_mixing(ky_mix,gamma_mix;kwargs...)
#
#    finds the maximum of gamma/ky spectrum vzf_out and kymax_out
#

    # uses variables from tjlf_dimensions, tglf_global, tglf_species
    kwargs=Dict(kwargs)

    #    initialize output of subroutine
    vzf_mix = 0.0
    kymax_mix = 0.0
    jmax_mix = 1
    xmin = 0.0

#
# find the local maximum of gamma_mix/ky_mix with the largest gamma_mix/ky_mix^2
#
    gammamax1= gamma_mix[1]
    kymax1 = ky_mix[1]
    testmax = 0.0

    ### kycut is defined with global variables in the python version, but rho_ion is a globalVar
    kycut=0.8/kwargs["rho_ion"]
    kymin = 0.0    

    
    kwargs["alpha_zf_in"] < 0.0 ? kymin = 0.173 * √(2.0) / kwargs["rho_ion"] :
    # saturation rules
    if sat_rule_in==2 || sat_rule_in==3
        kycut = kwargs["grad_r0_out"] * kycut
        kymin = kwargs["grad_r0_out"] * kymin
    end

    # find the low and high ky peaks of gamma/ky
    # nky = length(ky_mix) THIS LINE IS LEFTOVER FROM FORTRAN

    # go through ky_mix except for last value
    # if val less than kymin save the index
    # if next val is gte kymin and curr val is less than kycut save it as 1
    # update testmax and max index
    j1 = -1
    for j in range(1, length(ky_mix)-1)
        ### jmin not used
        ky_mix[j]<kymin ? jmin = j :
        
        if((ky_mix[j+1] >= kymin) && (ky_mix[j] <= kycut))
            j1=j
            kymax1 = ky_mix[j]
            testmax1 = gamma_mix[j]/kymax1
            # find maxes
            if(testmax1 > testmax)
                testmax = testmax1
                jmax_mix = j
            end
        end
    end

    ###### What do you do if jmax_mix = NOTHING
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
        if(f0 > f1)
# if f0>f1 then f1 is not a local maximum
            kymax1 = ky_mix[jmax1-1]
            gammamax1 = f0*kymax1
            if(kymax1 < kymin)
#interpolate to find the value of gammamax1 at kymin
                kymax1 = kymin
                xmin = (kymin - ky_mix[jmax1-1])/deltaky
                gammamax1 = (a + b*xmin + c*xmin*xmin)*kymin
            end
        end

        if(f0<f1)
#                if f0<f1 then f1>f2 due to the maximum search
# use the quadratic fit to refine the local maximum:
            xmax = -b/(2.0*c)
            xmin = 0.0
            if(ky_mix[jmax1-1]<kymin)
                xmin = (kymin - ky_mix[jmax1-1])/deltaky
            end

            if(xmax >= 1.0)
# if xmax >= 1    use f2 as the maximum
                kymax1 = ky_mix[jmax1+1]
                gammamax1 = f2*kymax1
            elseif(xmax <= xmin)    
                if(xmin > 0.0)
                    kymax1 = kymin
# use the quadratic fit to determine gammamax1 at kymin
                    gammamax1 = (a + b*xmin + c*xmin*xmin)*kymin
                elseif(xmax < 0.0)
# if xmax<=0 use f0 as the maximum
                    kymax1 = ky_mix[jmax1-1]
                    gammamax1 = f0*kymax1
                end
            else
# the conditions f0<f1<f2 and xmin<xmax<1 are satisfied
# use the quadratic fit to determine gammamax1 and kymax1
                kymax1 = ky_mix[jmax1-1]+deltaky*xmax
                gammamax1 = (a+b*xmax+c*xmax*xmax)*kymax1
            end #xmax tests
        end #f0 < f1
    end    # jmax_mix > 1
    vzf_mix = gammamax1/kymax1
    kymax_mix = kymax1
#            jmax_mix = jmax1
#            write(*,*)"get_zonal_mxing"
#            write(*,*)"xmax = ",xmax, "    xmin = ",xmin
#            write(*,*)"gammamax1 = ",gammamax1," kymax1 = ",kymax1, "    kymin = ",kymin
 
end
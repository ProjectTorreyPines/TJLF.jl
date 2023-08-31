

include("tjlf_modules.jl")

function get_zonal_mixing(nmix,ky_mix,gamma_mix,vzf_mix,kymax_mix,jmax_mix)
#
#  finds the maximum of gamma/ky spectrum vzf_out and kymax_out
#

  # uses variables from tjlf_dimensions, tglf_global, tglf_species

  use_kymin = false

  # INTEGER :: nmix,i,k,j,j1,j2,jmax1,jmin
  # REAL :: test,testmax,peakmax
  # REAL :: kx, kyhigh, kycut, kymin
  # REAL :: gammamax1,kymax1,testmax1,ky0,ky1,ky2
  # REAL :: f0,f1,f2,a,b,c,x1,deltaky,xmax,xmin
  # REAL :: vzf1, vzf2
  # REAL :: kymax_mix, vzf_mix
  # INTEGER :: jmax_mix, down
  # REAL, DIMENSION(nkym) :: gamma_mix, ky_mix

  #  initialize output of subroutine
  vzf_mix = 0.0
  kymax_mix = 0.0
  jmax_mix = 1
  xmin = 0.0

  tjlf_global.alpha_zf_in < 0.0 ? use_kymin = true :

#
# find the local maximum of gamma_mix/ky_mix with the largest gamma_mix/ky_mix^2
#
  gammamax1= gamma_mix[1]
  kymax1 = ky_mix[1]
  testmax = 0.0
  peakmax=0.0
  jmax_mix=1
  kycut=0.8/rho_ion
  kymin = 0.0
  jmin = 0
  
  use_kymin ? kymin = 0.173*√(2.0)/rho_ion :
  # saturation rules
  if sat_rule_in==2 || sat_rule_in==3
    kycut = grad_r0_out*kycut
    kymin = grad_r0_out*kymin
  end

# write(*,*)" kycut = ",kycut," kymin = ",kymin
  # find the low and high ky peaks of gamma/ky
  # nky = length(ky_mix) THIS LINE IS LEFTOVER FROM FORTRAN

  # go through ky_mix except for last value
  # if val less than kymin save the index
  # if next val is gte kymin and curr val is less than kycut save it as 1
  # update testmax and max index
  for j in range(1, length(ky_mix)-1)
    ### This seems redundant with the next if statement
    ky_mix[j]<kymin ? jmin = j :
    
    if((ky_mix[j+1] >= kymin) && (ky_mix[j] <= kycut))
      j1=j
      # redundant line because of line 81
      # kymax1 = ky_mix[j]
      testmax1 = gamma_mix[j]/kymax1
#     write(*,*)"j=",j,"ky = ",ky0," gamma_net = ",gamma_net(j)
      # find maxes
      if(testmax1 > testmax)
        testmax = testmax1
        jmax_mix = j
      end
    end
  end
  # if testmax is not updated a single time
  testmax==0.0 ? jmax_mix=j1 :
  # no unstable modes in range set kymax index to end of range
  # this is cut of at j1 since a maximum may not exist in the low-k range
  ### I don't like these line
  kymax1 = ky_mix[jmax_mix]
  gammamax1 = gamma_mix[jmax_mix]


  if(kymax1<kymin)
    kymax1 = kymin
# interpolate to find a more accurate low-k maximum gamma/ky
    gammamax1 = gamma_mix[1]
                +(gamma_mix[2]-gamma_mix[1])*(kymin-ky_mix[1])/(ky_mix[2]-ky_mix[1])
  end


#        write(*,*)" jmax_mix = ",jmax_mix,"  gammamax1 = ",gammamax1," kymax1 = ",kymax1
  if(jmax_mix>1 && jmax_mix < j1)
#        write(*,*)"refining low-k maximum"
# determine kymax1 and gammamax1 bounded by the tree points f0,f1,f2
# use a quadratic fit: f = a + b x + c x^2  to f = gamma/ky centered at jmax1
    jmax1 = jmax_mix
    f0 =  gamma_mix[jmax1-1]/ky_mix[jmax1-1]
    f1 =  gamma_mix[jmax1]/ky_mix[jmax1]
    f2 =  gamma_mix[jmax1+1]/ky_mix[jmax1+1]
    deltaky = ky_mix[jmax1+1]-ky_mix[jmax1-1]
    x1 = (ky_mix[jmax1]-ky_mix[jmax1-1])/deltaky
    a = f0
    b = (f1 - f0*(1-x1*x1)-f2*x1*x1)/(x1-x1*x1)
    c = f2 - f0 - b
#         write(*,*)"f0 = ",f0,"  f1 = ",f1,"  f2 = ",f2,"  x1 = ",x1
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
#        if f0<f1 then f1>f2 due to the maximum search
# use the quadratic fit to refine the local maximum:
      xmax = -b/(2.0*c)
      xmin = 0.0
      if(ky_mix[jmax1-1]<kymin)
        xmin = (kymin - ky_mix[jmax1-1])/deltaky
      end

      if(xmax >= 1.0)
# if xmax >= 1  use f2 as the maximum
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
  end  # jmax_mix > 1
  vzf_mix = gammamax1/kymax1
  kymax_mix = kymax1
#      jmax_mix = jmax1
#      write(*,*)"get_zonal_mxing"
#      write(*,*)"xmax = ",xmax, "  xmin = ",xmin
#      write(*,*)"gammamax1 = ",gammamax1," kymax1 = ",kymax1, "  kymin = ",kymin
 
end
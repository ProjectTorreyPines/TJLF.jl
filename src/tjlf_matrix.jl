using LinearAlgebra
import LinearAlgebra.LAPACK.syev!

include("tjlf_modules.jl")
include("tjlf_finiteLarmorRadius.jl")


function get_matrix(inputs::InputTJLF, ky)

    nbasis = inputs.NBASIS_MAX
    ns = inputs.NS

    ave = Ave{Float64}(ns, nbasis)
    aveH = AveH{Float64}(ns, nbasis)
    aveW = AveW{Float64}(ns, nbasis)
    aveK = AveK(ns,nbasis)

    ####### need gauher and gauss-hermite() functions for wx and h matrix
    FLR_xgrid!()
    get_ave!()

    if(inputs.VPAR_MODEL==0 && inputs.USE_BPER)
        ave.bpinv = inv(ave.bp)
        for i = 1:nbasis
            for j = 1:nbasis
                ave.p0inv[i,j] = 0.0
                ave.b0inv[i,j] = 0.0
                for k = 1:nbasis
                    ave.p0inv[i,j] = ave.p0inv[i,j]+ave.bpinv[i,k]*ave.b0[k,j]
                    ave.b0inv[i,j] = ave.b0inv[i,j]+ave.bpinv[i,k]*ave.p0[k,j]
                end
            end
        end
    else
        ave.p0inv = inv(ave.p0)
        if(inputs.USE_BPER) ave.b0inv = inv(ave.b0) end
    end

    modwd!()
    modkpar!()

    h_ratios!()

    if(inputs.LINSKER_FACTOR!=0.0) error("Haven't implemented :)") grad_ave_h() end

    ave_hp0!()
    if(inputs.BETAE>0.0)
        ave_hb0!()
        if(inputs.VPAR_MODEL==0) ave_hbp!() end
    end

    wd_h!()
    kpar_h!()
       
    if(inputs.GRADB_FACTOR!=0.0) error("Haven't implemented :)") gradB_h() end

end




#*************************************************************
#  begin computation of the FLR integrals
#*************************************************************
#  compute the FLR integrals at the hermite nodes
function FLR_xgrid!(ky::T, inputs::InputTJLF, ave::Ave, outputGeo::OutputGeometry, wx, h) where T<:Real

    nx = 2*inputs.NXGRID - 1
    ns = inputs.NS
    nbasis = inputs.NBASIS_MAX
    ns0 = 1
    if(inputs.ADIABATIC_ELEC) ns0 = 2 end

    b0x = outputGeo.b0x
    b2x = outputGeo.b2x

    fth=1.0
    hxn = Matrix{Float64}(undef, ns, nx) ### might not be floats
    for i = 1:nx
        for is =ns0:ns
            taus = inputs.SPECIES[is].TAUS
            mass = inputs.SPECIES[is].MASS
            zs = inputs.SPECIES[is].ZS

            bb = taus*mass*(ky/zs)^2 
            b = bb*b0x[i]/b2x[i]
            hn    = FLR_Hn(fth,b)
            hp1   = hn
            hp3   = FLR_dHp3(fth,b)+hn
            hr11  = 3.0*hp1
            hr13  = FLR_dHr13(fth,b)+(5.0/3.0)*hp1
            hr33  = FLR_dHr33(fth,b)+(5.0/3.0)*hp3
            hw113 = FLR_dHw113(fth,b)+(7.0/3.0)*hr11
            hw133 = FLR_dHw133(fth,b)+(7.0/3.0)*hr13
            hw333 = FLR_dHw333(fth,b)+(7.0/3.0)*hr33

            #***************************************************************
            #   compute the average h-bessel functions
            #***************************************************************
            for x = 1:nbasis
                for y = x:nbasis
                    ww = wx[i]*h[x,i]*h[y,i]
                    ave.hn[is,x,y]   = ave.hn[is,x,y]  + hn*ww
                    ave.hp1[is,x,y]  = ave.hp1[is,x,y]  + hp1*ww
                    ave.hp3[is,x,y]  = ave.hp3[is,x,y]  + hp3*ww
                    ave.hr11[is,x,y] = ave.hr11[is,x,y] + hr11*ww
                    ave.hr13[is,x,y] = ave.hr13[is,x,y] + hr13*ww
                    ave.hr33[is,x,y] = ave.hr33[is,x,y] + hr33*ww
                    ave.hw113[is,x,y] = ave.hw113[is,x,y] + hw113*ww
                    ave.hw133[is,x,y] = ave.hw133[is,x,y] + hw133*ww
                    ave.hw333[is,x,y] = ave.hw333[is,x,y] + hw333*ww

                    if(abs(ave.hn[is,x,y]).lt.zero_cut) ave.hn[is,x,y]=0.0 end
                    if(abs(ave.hp1[is,x,y]).lt.zero_cut) ave.hp1[is,x,y]=0.0 end
                    if(abs(ave.hp3[is,x,y]).lt.zero_cut) ave.hp3[is,x,y]=0.0 end
                    if(abs(ave.hr11[is,x,y]).lt.zero_cut) ave.hr11[is,x,y]=0.0 end
                    if(abs(ave.hr13[is,x,y]).lt.zero_cut) ave.hr13[is,x,y]=0.0 end
                    if(abs(ave.hr33[is,x,y]).lt.zero_cut) ave.hr33[is,x,y]=0.0 end
                    if(abs(ave.hw113[is,x,y]).lt.zero_cut) ave.hw113[is,x,y]=0.0 end
                    if(abs(ave.hw133[is,x,y]).lt.zero_cut) ave.hw133[is,x,y]=0.0 end
                    if(abs(ave.hw333[is,x,y]).lt.zero_cut) ave.hw333[is,x,y]=0.0 end

                    ave.hn[is,y,x]   = ave.hn[is,x,y] 
                    ave.hp1[is,j,i]  = ave.hp1[is,x,y]
                    ave.hp3[is,y,x]  = ave.hp3[is,x,y] 
                    ave.hr11[is,y,x] = ave.hr11[is,x,y] 
                    ave.hr13[is,y,x] = ave.hr13[is,x,y]
                    ave.hr33[is,y,x] = ave.hr33[is,x,y]
                    ave.hw113[is,y,x] = ave.hw113[is,x,y] 
                    ave.hw133[is,y,x] = ave.hw133[is,x,y]
                    ave.hw333[is,y,x] = ave.hw333[is,x,y]
                end
            end


        end
    end

    


end




#***********************************************************
#  compute  k-independent hermite basis averages
#***********************************************************
function get_ave!(ky::T,inputs::InputTJLF,ave::Ave,outputGeo::OutputGeometry,fts::Vector{T}, wx, h) where T<:Real

    zero_cut = 1.e-12
    ft2 = fts[1]^2

    vpar_model_in = inputs.VPAR_MODEL
    alpha_mach_in = inputs.ALPHA_MACH
    use_bper_in = inputs.USER_BPER

    betae_s = inputs.BETAE ##### not true for 'GENE' units
    debye_s = inputs.DEBYE #### different if UNITS is GENE
    debye_factor_in = inputs.DEBYE_FACTOR
    nx = 2*inputs.NXGRID - 1
    nbasis = inputs.NBASIS_MAX ##### might be wrong
    ns = inputs.NS
    ns0 = 1
    if(inputs.ADIABATIC_ELEC) ns0 = 2 end

    b0x = outputGeo.b0x
    wdx = outputGeo.wdx

    # fill the Poisson equation phi multiplier px0 x-grid array
    # note that pol is set in get_species
    pol = 0.0 ##### defined in  tglf_startup.f90
    for is = 1:ns
        taus = inputs.SPECIES[is].TAUS
        zs = inputs.SPECIES[is].ZS
        as = inputs.SPECIES[is].AS
        pol = pol +  zs^2*as/taus
    end
    U0 = 0.0 ### defined in startup.f90
    for is = 1:ns
        vpar_s = inputs.ALPHA_MACH*inputs.SIGN_IT*inputs.SPECIES[is].VPAR
        taus = inputs.SPECIES[is].TAUS
        as = inputs.SPECIES[is].AS
        zs = inputs.SPECIES[is].ZS
        
        U0 = U0 + as*vpar_s*zs^2/taus
    end


    p0x = Vector{Float64}(undef, nx)
    for i = 1:nx
        debye = debye_factor_in*b0x[i]*(ky*debye_s)^2
        p0x[i] = debye + pol
    end


    #  compute the guass-hermite intregrals
    for i = 1:nbasis        
        if(i<nbasis)
            ave.kpar[i,i+1] = √(i/2.0)
            ave.kpar[i+1,i] = -ave.kpar[i,i+1]
        end
       
        for j = i:nbasis
        
            for k = 1:nx
                ww = wx[k]*h[i,k]*h[j,k]     
                ave.wdh[i,j] = ave.wdh[i,j] + ww*wdx[k]
                ave.wdg[i,j] = ave.wdg[i,j] + ww*(wdx[k]+wdpx[k]*(1.0-ft2)/(1.0+ft2))
                ave.b0[i,j]  = ave.b0[i,j]  + ww*b0x[k]
                ave.p0inv[i,j] = ave.p0inv[i,j] + ww/p0x[k]
                ave.p0[i,j] = ave.p0[i,j] + ww*p0x[k]
                ave.kx[i,j] = ave.kx[i,j] + ww*kxx[k]
                ave.c_tor_par[i,j] = ave.c_tor_par[i,j] + ww*cx_tor_par[k]
                ave.c_tor_per[i,j] = ave.c_tor_per[i,j] + ww*cx_tor_per[k]
                ave.c_par_par[i,j] = ave.c_par_par[i,j] + ww*cx_par_par[k]
            end

            if(abs(ave.wdh[i,j])<zero_cut) ave.wdh[i,j] = 0.0 end
            if(abs(ave.wdg[i,j])<zero_cut) ave.wdg[i,j] = 0.0 end
            if(abs(ave.b0[i,j])<zero_cut) ave.b0[i,j] = 0.0 end
            if(abs(ave.lnB[i,j])<zero_cut) ave.lnB[i,j] = 0.0 end
            if(abs(ave.p0inv[i,j])<zero_cut) ave.p0inv[i,j] = 0.0 end
            if(abs(ave.p0[i,j])<zero_cut) ave.p0[i,j] = 0.0 end
            if(abs(ave.kx[i,j])<zero_cut) ave.kx[i,j] = 0.0 end
            if(abs(ave.c_tor_par[i,j])<zero_cut) ave.c_tor_par[i,j] = 0.0 end
            if(abs(ave.c_tor_per[i,j])<zero_cut) ave.c_tor_per[i,j] = 0.0 end
            if(abs(ave.c_par_par[i,j])<zero_cut) ave.c_par_par[i,j] = 0.0 end

            ave.wdh[j,i]       = ave.wdh[i,j]
            ave.wdg[j,i]       = ave.wdg[i,j]
            ave.b0[j,i]        = ave.b0[i,j]
            ave.lnB[j,i]       = ave.lnB[i,j]
            ave.p0inv[j,i]     = ave.p0inv[i,j]
            ave.p0[j,i]        = ave.p0[i,j]
            ave.kx[j,i]        = ave.kx[i,j]
            ave.c_tor_par[j,i] = ave.c_tor_par[i,j]
            ave.c_tor_per[j,i] = ave.c_tor_per[i,j]
            ave.c_par_par[j,i] = ave.c_par_par[i,j]
        end
    end

    for i = 1:nbasis
        for j = 1:nbasis
            gradB = 0.0
            for k = 1:nbasis
                gradB = gradB + ave.kpar[i,k]*ave.lnB[k,j] - ave.lnB[i,k]*ave.kpar[k,j]
            end
        ave.gradB[i,j] = gradB
        end
    end

    for is = ns0:ns
        if(nbasis==1)
            ave.kpar_eff[is,1,1] = -im/√(2)
            if(vpar_model_in==1)
                error("NOT IMPLEMENTED YET -DSUN")
                vpar_s = inputs.ALPHA_MACH*inputs.SIGN_IT*inputs.SPECIES[is].VPAR
                ave.kpar_eff[is,1,1] = ave.kpar_eff[is,1,1] + im*abs(alpha_mach_in*vpar_s)*ky*R_unit*q_unit*width_in*mass(is)/zs(is)
            end
        else
            for i = 1:nbasis
                for j = 1:nbasis
                    ave.kpar_eff[is,i,j] = ave.kpar[i,j]     
                    if(vpar_model_in==1 && i==j)
                        error("too lazy rn -DSUN")
                        vpar_s = inputs.ALPHA_MACH*inputs.SIGN_IT*inputs.SPECIES[is].VPAR
                        ave.kpar_eff[is,i,j] = ave.kpar_eff[is,i,j] - im*alpha_mach_in*ky*R_unit*q_unit*width_in*mass(is)/zs(is)
                    end
                end
            end
        end
    end

    if(vpar_model_in==0 && use_bper_in)
        betaU = 0.5*betae_s/(ky*ky)
        betaU = betaU*U0^2
        for i = 1:nbasis
            for j = 1:nbasis
                ave.bp[i,j] = 0.0
                for k = 1:nbasis
                    ave.bp[i,j] = ave.bp[i,j] + ave.b0[i,k]*ave.p0[k,j]
                end
                if(i==j) ave.bp[i,j] = ave.bp[i,j] + betaU end
            end
        end
    end

end



#***************************************************************
#   compute the matricies modwdh and modwdg
#***************************************************************
function modwd!(inputs::InputTJLF,ave::Ave)

    nbasis = inputs.NBASIS_MAX
    wd_zero_in = inputs.WD_ZERO
    nm = nbasis

    if(nm==1)
        ave.modwdh[1,1] = abs(ave.wdh[1,1])
        ave.modwdg[1,1] = abs(ave.wdg[1,1])
        return
    end

    # find the eigenvalues of ave.wdh
    a = deepcopy(ave.wdh)
    w,a = syev!('V','U',a)
    for k = 1:nm
        if(abs(w[k])<wd_zero_in)
            if(w[k]>=0.0)
                w[k] = wd_zero_in
            else
                w[k] = -wd_zero_in
            end
        end
    end
    # compute ave.modwd and recompute ave.wd with regularized eigenvalues
    # note that the DSYEV normalized eigenvectors are now in a[i,j]
    w = Diagonal(w)
    ave.modwdh = a * abs.(w) * a
    ave.wdh = a * w * a ### excessive

    # find the eigenvalues of ave.wdg
    a = deepcopy(ave.wdg)
    w,a = syev!('V','U',a)
    for k = 1:nm
        if(abs(w[k])<wd_zero_in)
            if(w[k]>=0.0)
                w[k] = wd_zero_in
            else
                w[k] = -wd_zero_in
            end
        end
    end
    # compute ave.modwd and recompute ave.wd with regularized eigenvalues
    # note that the DSYEV normalized eigenvectors are now in a[i,j]
    w = Diagonal(w)
    ave.modwdg = a * abs.(w) * a
    ave.wdg = a * w * a ### excessive

end




#***************************************************************
#   compute the matrix mod_kpar
#***************************************************************
function modkpar!(inputs::InputTJLF,ave::Ave)

    vpar_model_in = inputs.VPAR_MODEL
    ns0 = 1
    if(inputs.ADIABATIC_ELEC) ns0 = 2 end
    nbasis = inputs.NBASIS_MAX
    nm = nbasis

    if(nm==1)
        for is = ns0:ns
            ave.modkpar_eff[is,1,1] = im*ave.kpar_eff[is,1,1]
        end
    else

    # find the eigenvalues and eigenvectors
        if(vpar_model_in!=1)   
            a = im*deepcopy(ave.kpar)

            w,a = syev!('V','U',a)
            w = Diagonal(w)
            b = a * abs.(w) * conj.(a)
            ave.modkpar = real.(b)
            for is = ns0:ns
                ave.modkpar_eff[is,:,:] .= b
            end

        else 
            for is = ns0:ns
                a = im*deepcopy(ave.kpar_eff[is,:,:])
                w,a = syev!('V','U',a)
                b = a * abs.(w) * conj.(a)
                
                ave.modkpar_eff[is,:,:] .= b
            end 
        end

    end
end





#***************************************************************
#   compute the h_ratios
#***************************************************************
function h_ratios!(inputs::InputTJLF, aveH::AveH, )

    ns = inputs.NS
    ns0 = 1
    if(inputs.ADIABATIC_ELEC) ns0 = 2 end
    hninv = inv(aveH.hn)
    hp1inv = inv(aveH.hp1)
    hp3inv = inv(aveH.hp3)
    # compute matrix ratios
    for is = ns0:ns
        aveH.ht1[is,:,:] .= aveH.hp1 * hninv
        aveH.ht3[is,:,:] .= aveH.hp3 * hninv
        aveH.hu1[is,:,:] .= aveH.hr11 * hp1inv
        aveH.hu3[is,:,:] .= aveH.hr13 * hp1inv
        aveH.hu33[is,:,:] .= aveH.hr33 * hp3inv

        aveH.hu3ht1[is,:,:] .= aveH.hu3[is,:,:] * aveH.ht1[is,:,:]
        aveH.hu33ht1[is,:,:] .= aveH.hu33[is,:,:] * aveH.ht1[is,:,:]
        aveH.hu3ht3[is,:,:] .= aveH.hu3[is,:,:] * aveH.ht3[is,:,:]
        aveH.hu33ht3[is,:,:] .= aveH.hu33[is,:,:] * aveH.ht3[is,:,:]
    end

end



#***************************************************************
#   compute the products ave_h*ave_p0inv
#***************************************************************
function ave_hp0!(inputs::InputTJLF,aveH::AveH)

    ns = inputs.NS
    ns0 = 1
    if(inputs.ADIABATIC_ELEC) ns0 = 2 end
    # compute matrix products
    for is = ns0:ns
        aveH.hnp0[is,:,:]    .= aveH.hn[is,:,:] * aveH.p0inv
        aveH.hp1p0[is,:,:]   .= aveH.hp1[is,:,:]*aveH.p0inv
        aveH.hp3p0[is,:,:]   .= aveH.hp3[is,:,:]*aveH.p0inv
        aveH.hr11p0[is,:,:]  .= aveH.hr11[is,:,:]*aveH.p0inv
        aveH.hr13p0[is,:,:]  .= aveH.hr13[is,:,:]*aveH.p0inv
        aveH.hr33p0[is,:,:]  .= aveH.hr33[is,:,:]*aveH.p0inv     
        
        aveH.c_tor_par_hp1p0[is,:,:]  .= aveH.c_tor_par * aveH.hp1p0[is,:,:]
        aveH.c_tor_par_hr11p0[is,:,:] .= aveH.c_tor_par * aveH.hr11p0[is,:,:]
        aveH.c_tor_par_hr13p0[is,:,:] .= aveH.c_tor_par * aveH.hr13p0[is,:,:]
    end

end


#***************************************************************
#   compute the products ave_h*ave_b0inv
#***************************************************************
function ave_hb0!(inputs::InputTJLF,ave::Ave,aveH::AveH)

    ns = inputs.NS
    ns0 = 1
    if(inputs.ADIABATIC_ELEC) ns0 = 2 end
    # compute matrix products
    for is = ns0:ns
        aveH.hnb0[is,:,:]     .= aveH.hn[is,:,:] * ave.b0inv
        aveH.hp1b0[is,:,:]    .= aveH.hp1[is,:,:] * ave.b0inv
        aveH.hp3b0[is,:,:]    .= aveH.hp3[is,:,:] * ave.b0inv
        aveH.hr11b0[is,:,:]   .= aveH.hr11[is,:,:] * ave.b0inv
        aveH.hr13b0[is,:,:]   .= aveH.hr13[is,:,:] * ave.b0inv
        aveH.hr33b0[is,:,:]   .= aveH.hr33[is,:,:] * ave.b0inv
        aveH.hw113b0[is,:,:]  .= aveH.hw113[is,:,:] * ave.b0inv
        aveH.hw133b0[is,:,:]  .= aveH.hw133[is,:,:] * ave.b0inv
        aveH.hw333b0[is,:,:]  .= aveH.hw333[is,:,:] * ave.b0inv
    end
end



#***************************************************************
#   compute the products ave_h*ave_bpinv
#***************************************************************
function ave_hbp!(inputs::InputTJLF,ave::Ave,aveH::AveH)

    ns = inputs.NS
    ns0 = 1
    if(inputs.ADIABATIC_ELEC) ns0 = 2 end

    for is = ns0:ns
        aveH.hnbp[is,:,:]     .= aveH.hn[is,:,:] * ave.bpinv
        aveH.hp1bp[is,:,:]    .= aveH.hp1[is,:,:] * ave.bpinv
        aveH.hp3bp[is,:,:]    .= aveH.hp3[is,:,:] * ave.bpinv
        aveH.hr11bp[is,:,:]   .= aveH.hr11[is,:,:] * ave.bpinv
        aveH.hr13bp[is,:,:]   .= aveH.hr13[is,:,:] * ave.bpinv
        aveH.hr33bp[is,:,:]   .= aveH.hr33[is,:,:] * ave.bpinv
        aveH.hw113bp[is,:,:]  .= aveH.hw113[is,:,:] * ave.bpinv
        aveH.hw133bp[is,:,:]  .= aveH.hw133[is,:,:] * ave.bpinv
        aveH.hw333bp[is,:,:]  .= aveH.hw333[is,:,:] * ave.bpinv
    end       
end


#***************************************************************
#   compute the products ave_wd*ave_h
#***************************************************************
function wd_h!(inputs::InputTJLF,ave::Ave,aveH::AveH)

    use_bper_in = inputs.USE_BPER
    vpar_model_in = inputs.VPAR_MODEL
    ns = inputs.NS
    ns0 = 1
    if(inputs.ADIABATIC_ELEC) ns0 = 2 end

    for is = ns0:ns

        ave_wdhp1p0[is,:,:]  .= ave.wdh * aveH.hp1p0[is,:,:]
        ave_wdhr11p0[is,:,:] .= ave.wdh * aveH.hr11p0[is,:,:]
        ave_wdhr13p0[is,:,:] .= ave.wdh * aveH.hr13p0[is,:,:]
        ave_wdht1[is,:,:]    .= ave.wdh * aveH.ht1[is,:,:]
        ave_wdht3[is,:,:]    .= ave.wdh * aveH.ht3[is,:,:]
        ave_wdhu1[is,:,:]    .= ave.wdh * aveH.hu1[is,:,:]
        ave_wdhu3[is,:,:]    .= ave.wdh * aveH.hu3[is,:,:]
        ave_wdhu3ht1[is,:,:] .= ave.wdh * aveH.hu3ht1[is,:,:]
        ave_wdhu3ht3[is,:,:] .= ave.wdh * aveH.hu3ht3[is,:,:]
        ave_wdhu33[is,:,:]   .= ave.wdh * aveH.hu33[is,:,:]
        ave_wdhu33ht1[is,:,:].= ave.wdh * aveH.hu33ht1[is,:,:]
        ave_wdhu33ht3[is,:,:].= ave.wdh * aveH.hu33ht3[is,:,:]
        ave_modwdht1[is,:,:] .= ave.modwdh * aveH.ht1[is,:,:]
        ave_modwdht3[is,:,:] .= ave.modwdh * aveH.ht3[is,:,:]
        ave_modwdhu1[is,:,:] .= ave.modwdh * aveH.hu1[is,:,:]
        ave_modwdhu3[is,:,:] .= ave.modwdh * aveH.hu3[is,:,:]
        ave_modwdhu3ht1[is,:,:] .= ave.modwdh * aveH.hu3ht1[is,:,:]
        ave_modwdhu3ht3[is,:,:] .= ave.modwdh * aveH.hu3ht3[is,:,:]
        ave_modwdhu33[is,:,:]   .= ave.modwdh * aveH.hu33[is,:,:]
        ave_modwdhu33ht1[is,:,:].= ave.modwdh * aveH.hu33ht1[is,:,:]
        ave_modwdhu33ht3[is,:,:].= ave.modwdh * aveH.hu33ht3[is,:,:]
        
        if(use_bper_in)
            ave_wdhp1b0[is,:,:]  .= ave.wdh * aveH.hp1b0[is,:,:]
            ave_wdhr11b0[is,:,:] .= ave.wdh * aveH.hr11b0[is,:,:]
            ave_wdhr13b0[is,:,:] .= ave.wdh * aveH.hr13b0[is,:,:]
            if(vpar_model_in==0)
                ave_wdhp1bp[is,:,:]  .= ave.wdh * aveH.hp1bp[is,:,:]
                ave_wdhr11bp[is,:,:] .= ave.wdh * aveH.hr11bp[is,:,:]
                ave_wdhr13bp[is,:,:] .= ave.wdh * aveH.hr13bp[is,:,:]
            end
        end
    end
end





#***************************************************************
#   compute the products ave_kpar*ave_h and ave_modkpar*ave_h
#***************************************************************
function kpar_h!(inputs::InputTJLF,ave::Ave,aveH::AveH,aveK::AveK)

    use_bper_in = inputs.USE_BPER
    vpar_model_in = inputs.VPAR_MODEL
    ns = inputs.NS
    ns0 = 1
    if(inputs.ADIABATIC_ELEC) ns0 = 2 end

    for is = ns0:ns
        ave_kparhnp0[is,:,:]  .= aveH.hnp0[is,:,:] * ave.kpar
        ave_kparhp1p0[is,:,:] .= aveH.hp1p0[is,:,:] * ave.kpar
        ave_kparhp3p0[is,:,:] .= aveH.hp3p0[is,:,:] * ave.kpar
        ave_kparhu1[is,:,:]   .= aveH.hu1[is,:,:] * ave.kpar_eff[is,:,:]
        ave_kparhu3[is,:,:]   .= aveH.hu3[is,:,:] * ave.kpar_eff[is,:,:]
        ave_kparht1[is,:,:]   .= ave.kpar_eff[is,:,:] * aveH.ht1[is,:,:]
        ave_kparht3[is,:,:]   .= ave.kpar_eff[is,:,:] * aveH.ht3[is,:,:]
        ave_modkparhu1[is,:,:].= ave.modkpar_eff[is,:,:] * aveH.hu1[is,:,:]
        ave_modkparhu3[is,:,:].= ave.modkpar_eff[is,:,:] * aveH.hu3[is,:,:]
        if(use_bper_in)
            ave_kparhp1b0[is,:,:]  .= aveH.hp1b0[is,:,:] * ave.kpar
            ave_kparhr11b0[is,:,:] .= aveH.hr11b0[is,:,:] * ave.kpar
            ave_kparhr13b0[is,:,:] .= aveH.hr13b0[is,:,:] * ave.kpar
            if(vpar_model_in==0)
                ave_kparhnbp[is,:,:]  .= aveH.hnbp[is,:,:] * ave.kpar
                ave_kparhp3bp[is,:,:] .= aveH.hp3bp[is,:,:] * ave.kpar
                ave_kparhp1bp[is,:,:] .= aveH.hp1bp[is,:,:] * ave.kpar
                ave_kparhr11bp[is,:,:].= aveH.hr11bp[is,:,:] * ave.kpar
                ave_kparhr13bp[is,:,:].= aveH.hr13bp[is,:,:] * ave.kpar
            end
        end
    end
       
end
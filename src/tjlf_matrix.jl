using LinearAlgebra
import LinearAlgebra.LAPACK.syev!

include("tjlf_modules.jl")
include("tjlf_finiteLarmorRadius.jl")


function get_matrix(inputs::InputTJLF{T}, outputGeo::OutputGeometry{T}, outputHermite::OutputHermite{T}, 
                    ky::T, 
                    nbasis::Int) where T<:Real

    ns = inputs.NS

    ave = Ave{Float64}(ns, nbasis)
    aveH = AveH{Float64}(ns, nbasis)
    aveWH = AveWH{Float64}(ns, nbasis)
    aveKH = AveKH(ns, nbasis)

    aveG = AveG{Float64}(ns, nbasis)
    aveWG = AveWG{Float64}(ns, nbasis)
    aveKG = AveKG(ns, nbasis)

    FLR_xgrid!(inputs, outputGeo, outputHermite, aveH, aveG, ky, nbasis)
    get_ave!(inputs, outputGeo, outputHermite, ave, nbasis, ky)

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

    modwd!(inputs, ave)
    modkpar!(inputs, ave)

    h_ratios!(inputs, aveH)

    if(inputs.LINSKER_FACTOR!=0.0) error("Haven't implemented :)") end #grad_ave_h() end

    ave_hp0!(inputs, ave, aveH)
    if(inputs.BETAE>0.0)
        ave_hb0!(inputs, ave, aveH)
        if(inputs.VPAR_MODEL==0) ave_hbp!(inputs, ave, aveH) end
    end

    wd_h!(inputs, ave, aveH, aveWH)
    kpar_h!(inputs, ave, aveH, aveKH)
       
    if(inputs.GRADB_FACTOR!=0.0) error("Haven't implemented :)") end #gradB_h() end


    nroot = 15 ###### hardcoded
    if(nroot>6)
        g_ratios!(inputs, aveG)
        if(inputs.LINSKER_FACTOR!=0) error("Have't implemented :)") end #grad_ave_g
        ave_gp0!(inputs, ave, aveG)
        if(inputs.BETAE>0.0)
            ave_gb0!(inputs, ave, aveG)
            if(inputs.VPAR_MODEL==0) ave_gbp!(inputs, ave, aveG) end
        end
        wd_g!(inputs, ave, aveG, aveWG)
        kpar_g!(inputs, ave, aveG, aveKG)
        if(inputs.GRADB_FACTOR!=0.0) error("Haven't implemented :)") end #gradB_g
    end

    return ave, aveH, aveWH, aveKH, aveG, aveWG, aveKG
end




#*************************************************************
#  begin computation of the FLR integrals
#*************************************************************
#  compute the FLR integrals at the hermite nodes
function FLR_xgrid!(inputs::InputTJLF{T}, outputGeo::OutputGeometry{T}, outputHermite::OutputHermite{T}, 
    aveH::AveH{T}, aveG::AveG{T}, ky::T, nbasis::Int) where T<:Real

    zero_cut = 1.e-12
    nx = 2*inputs.NXGRID - 1
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    b0x = outputGeo.b0x
    b2x = outputGeo.b2x
    wx = outputHermite.wx
    h = outputHermite.h[1:nbasis,:]

    
    hxn = Matrix{Float64}(undef, ns, nx)
    hxp1 = Matrix{Float64}(undef, ns, nx)
    hxp3 = Matrix{Float64}(undef, ns, nx)
    hxr11 = Matrix{Float64}(undef, ns, nx)
    hxr13 = Matrix{Float64}(undef, ns, nx)
    hxr33 = Matrix{Float64}(undef, ns, nx)
    hxw113 = Matrix{Float64}(undef, ns, nx)
    hxw133 = Matrix{Float64}(undef, ns, nx)
    hxw333 = Matrix{Float64}(undef, ns, nx)
    nroot = 15 ###### hardcoded
    if(nroot>6)
        gxn = Matrix{Float64}(undef, ns, nx)
        gxp1 = Matrix{Float64}(undef, ns, nx)
        gxp3 = Matrix{Float64}(undef, ns, nx)
        gxr11 = Matrix{Float64}(undef, ns, nx)
        gxr13 = Matrix{Float64}(undef, ns, nx)
        gxr33 = Matrix{Float64}(undef, ns, nx)
        gxw113 = Matrix{Float64}(undef, ns, nx)
        gxw133 = Matrix{Float64}(undef, ns, nx)
        gxw333 = Matrix{Float64}(undef, ns, nx)
    end

    fth = 1.0
    for i = 1:nx
        for is = ns0:ns
            taus = inputs.TAUS[is]
            mass = inputs.MASS[is]
            zs = inputs.ZS[is]

            #***************************************************************
            #   compute the average h-bessel functions
            #***************************************************************
            bb = taus*mass*(ky/zs)^2 
            b = bb*b0x[i]/b2x[i]
            hxn[is,i]    = FLR_Hn(fth,b)
            hxp1[is,i]   = hxn[is,i]
            hxp3[is,i]   = FLR_dHp3(fth,b) + hxn[is,i]
            hxr11[is,i]  = 3.0 * hxp1[is,i]
            hxr13[is,i]  = FLR_dHr13(fth,b) + (5/3)*hxp1[is,i]
            hxr33[is,i]  = FLR_dHr33(fth,b) + (5/3)*hxp3[is,i]
            hxw113[is,i] = FLR_dHw113(fth,b)+ (7/3)*hxr11[is,i]
            hxw133[is,i] = FLR_dHw133(fth,b)+ (7/3)*hxr13[is,i]
            hxw333[is,i] = FLR_dHw333(fth,b)+ (7/3)*hxr33[is,i]
            
            #***************************************************************
            #   compute the average g- bessel functions
            #***************************************************************
            if(nroot>6)
                fts = outputGeo.fts
                ftx::Float64 = fts[is]
                ft2 = ftx^2
                gxn[is,i]    = FLR_Hn(ftx,b)
                gxp1[is,i]   = FLR_dHp1(ftx,b)  + ft2*gxn[is,i]
                gxp3[is,i]   = FLR_dHp3(ftx,b)  + gxn[is,i]
                gxr11[is,i]  = FLR_dHr11(ftx,b) + 3*ft2*gxp1[is,i]
                gxr13[is,i]  = FLR_dHr13(ftx,b) + (5/3)*gxp1[is,i]
                gxr33[is,i]  = FLR_dHr33(ftx,b) + (5/3)*gxp3[is,i]
                gxw113[is,i] = FLR_dHw113(ftx,b)+ (7/3)*gxr11[is,i]
                gxw133[is,i] = FLR_dHw133(ftx,b)+ (7/3)*gxr13[is,i]
                gxw333[is,i] = FLR_dHw333(ftx,b)+ (7/3)*gxr33[is,i]

            end
        end
    end

    for is = ns0:ns
        aveH.hn[is,:,:]    .= h[1:nbasis,:] * Diagonal(hxn[is,:]   .*wx)   * h'
        aveH.hp1[is,:,:]   .= h[1:nbasis,:] * Diagonal(hxp1[is,:]  .*wx)   * h'
        aveH.hp3[is,:,:]   .= h[1:nbasis,:] * Diagonal(hxp3[is,:]  .*wx)   * h'
        aveH.hr11[is,:,:]  .= h[1:nbasis,:] * Diagonal(hxr11[is,:] .*wx)   * h'
        aveH.hr13[is,:,:]  .= h[1:nbasis,:] * Diagonal(hxr13[is,:] .*wx)   * h'
        aveH.hr33[is,:,:]  .= h[1:nbasis,:] * Diagonal(hxr33[is,:] .*wx)   * h'
        aveH.hw113[is,:,:] .= h[1:nbasis,:] * Diagonal(hxw113[is,:].*wx)   * h'
        aveH.hw133[is,:,:] .= h[1:nbasis,:] * Diagonal(hxw133[is,:].*wx)   * h'
        aveH.hw333[is,:,:] .= h[1:nbasis,:] * Diagonal(hxw333[is,:].*wx)   * h'

        aveH.hn[abs.(aveH.hn) .< zero_cut] .= 0
        aveH.hp1[abs.(aveH.hp1) .< zero_cut] .= 0
        aveH.hp3[abs.(aveH.hp3) .< zero_cut] .= 0
        aveH.hr11[abs.(aveH.hr11) .< zero_cut] .= 0
        aveH.hr13[abs.(aveH.hr13) .< zero_cut] .= 0
        aveH.hr33[abs.(aveH.hr33) .< zero_cut] .= 0
        aveH.hw113[abs.(aveH.hw113) .< zero_cut] .= 0
        aveH.hw133[abs.(aveH.hw133) .< zero_cut] .= 0
        aveH.hw333[abs.(aveH.hw333) .< zero_cut] .= 0
        
        nroot = 15 ####### hardcoded
        if(nroot>6)
            aveG.gn[is,:,:]    .= h * Diagonal(gxn[is,:]   .*wx)   * h'
            aveG.gp1[is,:,:]   .= h * Diagonal(gxp1[is,:]  .*wx)   * h'
            aveG.gp3[is,:,:]   .= h * Diagonal(gxp3[is,:]  .*wx)   * h'
            aveG.gr11[is,:,:]  .= h * Diagonal(gxr11[is,:] .*wx)   * h'
            aveG.gr13[is,:,:]  .= h * Diagonal(gxr13[is,:] .*wx)   * h'
            aveG.gr33[is,:,:]  .= h * Diagonal(gxr33[is,:] .*wx)   * h'
            aveG.gw113[is,:,:] .= h * Diagonal(gxw113[is,:].*wx)   * h'
            aveG.gw133[is,:,:] .= h * Diagonal(gxw133[is,:].*wx)   * h'
            aveG.gw333[is,:,:] .= h * Diagonal(gxw333[is,:].*wx)   * h'
            
            aveG.gn[abs.(aveG.gn) .< zero_cut] .= 0
            aveG.gp1[abs.(aveG.gp1) .< zero_cut] .= 0
            aveG.gp3[abs.(aveG.gp3) .< zero_cut] .= 0
            aveG.gr11[abs.(aveG.gr11) .< zero_cut] .= 0
            aveG.gr13[abs.(aveG.gr13) .< zero_cut] .= 0
            aveG.gr33[abs.(aveG.gr33) .< zero_cut] .= 0
            aveG.gw113[abs.(aveG.gw113) .< zero_cut] .= 0
            aveG.gw133[abs.(aveG.gw133) .< zero_cut] .= 0
            aveG.gw333[abs.(aveG.gw333) .< zero_cut] .= 0
        end
    end
end




#***********************************************************
#  compute  k-independent hermite basis averages
#***********************************************************
function get_ave!(inputs::InputTJLF{T},outputGeo::OutputGeometry{T},outputHermite::OutputHermite{T},ave::Ave{T}, 
    nbasis::Int, ky::T) where T<:Real

    zero_cut = 1.e-12
    fts = outputGeo.fts 
    ft2 = fts[1]^2

    width_in = inputs.WIDTH

    vpar_model_in = inputs.VPAR_MODEL
    alpha_mach_in = inputs.ALPHA_MACH
    sign_it_in = inputs.SIGN_IT
    use_bper_in = inputs.USE_BPER

    betae_s = inputs.BETAE ##### not true for 'GENE' units
    debye_s = inputs.DEBYE #### different if UNITS is GENE
    debye_factor_in = inputs.DEBYE_FACTOR
    nx = 2*inputs.NXGRID - 1
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    taus = inputs.TAUS
    vpar = inputs.VPAR
    zs = inputs.ZS
    as = inputs.AS

    b0x = outputGeo.b0x
    wdx = outputGeo.wdx
    wdpx = outputGeo.wdpx
    kxx = outputGeo.kxx
    cx_tor_par = outputGeo.cx_tor_par
    cx_tor_per = outputGeo.cx_tor_per
    cx_par_par = outputGeo.cx_par_par
    wx = outputHermite.wx
    h = outputHermite.h[1:nbasis,:]

    # fill the Poisson equation phi multiplier px0 x-grid array
    # note that pol is set in get_species
    pol = 0.0 ##### defined in  tglf_startup.f90
    pol = sum(zs.^2 .* as./taus)
    U0 = sum((alpha_mach_in*sign_it_in).*vpar.*zs.^2 .* as./taus) ### defined in startup.f90


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

    ave.gradB .= (ave.kpar*ave.lnB) .- (ave.lnB*ave.kpar)

    for is = ns0:ns
        if(nbasis==1)
            ave.kpar_eff[is,1,1] = -im/√(2)
            if(vpar_model_in==1)
                error("NOT IMPLEMENTED YET -DSUN")
                vpar_s = inputs.ALPHA_MACH*inputs.SIGN_IT*inputs.VPAR[is]
                ave.kpar_eff[is,1,1] = ave.kpar_eff[is,1,1] + im*abs(alpha_mach_in*vpar_s)*ky*R_unit*q_unit*width_in*mass(is)/zs(is)
            end
        else
            for i = 1:nbasis
                for j = 1:nbasis
                    ave.kpar_eff[is,i,j] = ave.kpar[i,j]     
                    if(vpar_model_in==1 && i==j)
                        error("too lazy rn -DSUN")
                        vpar_s = inputs.ALPHA_MACH*inputs.SIGN_IT*inputs.VPAR[is]
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
function modwd!(inputs::InputTJLF{T},ave::Ave{T}) where T<:Real

    wd_zero_in = inputs.WD_ZERO

    if(size(ave.modwdh)[2]==1)
        ave.modwdh[1,1] = abs(ave.wdh[1,1])
        ave.modwdg[1,1] = abs(ave.wdg[1,1])
        return
    end

    # find the eigenvalues of ave.wdh
    a = deepcopy(ave.wdh)
    w,a = syev!('V','U',a)
    for k = 1:size(ave.modwdh)[2]
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
    ave.modwdh = a * abs.(w) * transpose(a)
    ave.wdh = a * w * transpose(a) ### excessive

    # find the eigenvalues of ave.wdg
    a = deepcopy(ave.wdg)
    w,a = syev!('V','U',a)
    for k = 1:size(ave.modwdh)[2]
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
    ave.modwdg = a * abs.(w) * transpose(a)
    ave.wdg = a * w * transpose(a) ### excessive

end




#***************************************************************
#   compute the matrix mod_kpar
#***************************************************************
function modkpar!(inputs::InputTJLF{T},ave::Ave{T}) where T<:Real

    vpar_model_in = inputs.VPAR_MODEL
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    if(size(ave.modkpar_eff)[2]==1)
        for is = ns0:ns
            ave.modkpar_eff[is,1,1] = im*ave.kpar_eff[is,1,1]
        end
    else
    # find the eigenvalues and eigenvectors
        if(vpar_model_in!=1)   
            a = im*deepcopy(ave.kpar)

            w,a = syev!('V','U',a)
            w = Diagonal(w)
            b = a * abs.(w) * transpose(conj.(a))
            ave.modkpar = real.(b)
            for is = ns0:ns
                ave.modkpar_eff[is,:,:] .= b
            end

        else 
            for is = ns0:ns
                a = im*deepcopy(ave.kpar_eff[is,:,:])
                w,a = syev!('V','U',a)
                b = a * abs.(w) * transpose(conj.(a))
                
                ave.modkpar_eff[is,:,:] .= b
            end 
        end

    end
end





#***************************************************************
#   compute the h_ratios
#***************************************************************
function h_ratios!(inputs::InputTJLF{T}, aveH::AveH{T}) where T<:Real

    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    # compute matrix ratios
    for is = ns0:ns
        hninv = inv(aveH.hn[is,:,:])
        hp1inv = inv(aveH.hp1[is,:,:])
        hp3inv = inv(aveH.hp3[is,:,:])

        aveH.ht1[is,:,:] .= aveH.hp1[is,:,:] * hninv
        aveH.ht3[is,:,:] .= aveH.hp3[is,:,:] * hninv
        aveH.hu1[is,:,:] .= aveH.hr11[is,:,:] * hp1inv
        aveH.hu3[is,:,:] .= aveH.hr13[is,:,:] * hp1inv
        aveH.hu33[is,:,:] .= aveH.hr33[is,:,:] * hp3inv

        aveH.hu3ht1[is,:,:]  .= aveH.hu3[is,:,:]  * aveH.ht1[is,:,:]
        aveH.hu33ht1[is,:,:] .= aveH.hu33[is,:,:] * aveH.ht1[is,:,:]
        aveH.hu3ht3[is,:,:]  .= aveH.hu3[is,:,:]  * aveH.ht3[is,:,:]
        aveH.hu33ht3[is,:,:] .= aveH.hu33[is,:,:] * aveH.ht3[is,:,:]
    end

end



#***************************************************************
#   compute the products ave_h*ave_p0inv
#***************************************************************
function ave_hp0!(inputs::InputTJLF{T},ave::Ave{T},aveH::AveH{T}) where T<:Real

    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    # compute matrix products
    for is = ns0:ns
        aveH.hnp0[is,:,:]    .= aveH.hn[is,:,:]*ave.p0inv
        aveH.hp1p0[is,:,:]   .= aveH.hp1[is,:,:]*ave.p0inv
        aveH.hp3p0[is,:,:]   .= aveH.hp3[is,:,:]*ave.p0inv
        aveH.hr11p0[is,:,:]  .= aveH.hr11[is,:,:]*ave.p0inv
        aveH.hr13p0[is,:,:]  .= aveH.hr13[is,:,:]*ave.p0inv
        aveH.hr33p0[is,:,:]  .= aveH.hr33[is,:,:]*ave.p0inv     
        
        aveH.c_tor_par_hp1p0[is,:,:]  .= ave.c_tor_par * aveH.hp1p0[is,:,:]
        aveH.c_tor_par_hr11p0[is,:,:] .= ave.c_tor_par * aveH.hr11p0[is,:,:]
        aveH.c_tor_par_hr13p0[is,:,:] .= ave.c_tor_par * aveH.hr13p0[is,:,:]
    end

end


#***************************************************************
#   compute the products ave_h*ave_b0inv
#***************************************************************
function ave_hb0!(inputs::InputTJLF{T},ave::Ave{T},aveH::AveH{T}) where T<:Real

    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
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
function ave_hbp!(inputs::InputTJLF{T},ave::Ave{T},aveH::AveH{T}) where T<:Real

    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

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
function wd_h!(inputs::InputTJLF{T},ave::Ave{T},aveH::AveH{T},aveWH::AveWH{T}) where T<:Real

    use_bper_in = inputs.USE_BPER
    vpar_model_in = inputs.VPAR_MODEL
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    for is = ns0:ns

        aveWH.wdhp1p0[is,:,:]  .= ave.wdh * aveH.hp1p0[is,:,:]
        aveWH.wdhr11p0[is,:,:] .= ave.wdh * aveH.hr11p0[is,:,:]
        aveWH.wdhr13p0[is,:,:] .= ave.wdh * aveH.hr13p0[is,:,:]
        aveWH.wdht1[is,:,:]    .= ave.wdh * aveH.ht1[is,:,:]
        aveWH.wdht3[is,:,:]    .= ave.wdh * aveH.ht3[is,:,:]
        aveWH.wdhu1[is,:,:]    .= ave.wdh * aveH.hu1[is,:,:]
        aveWH.wdhu3[is,:,:]    .= ave.wdh * aveH.hu3[is,:,:]
        aveWH.wdhu3ht1[is,:,:] .= ave.wdh * aveH.hu3ht1[is,:,:]
        aveWH.wdhu3ht3[is,:,:] .= ave.wdh * aveH.hu3ht3[is,:,:]
        aveWH.wdhu33[is,:,:]   .= ave.wdh * aveH.hu33[is,:,:]
        aveWH.wdhu33ht1[is,:,:].= ave.wdh * aveH.hu33ht1[is,:,:]
        aveWH.wdhu33ht3[is,:,:].= ave.wdh * aveH.hu33ht3[is,:,:]
        aveWH.modwdht1[is,:,:] .= ave.modwdh * aveH.ht1[is,:,:]
        aveWH.modwdht3[is,:,:] .= ave.modwdh * aveH.ht3[is,:,:]
        aveWH.modwdhu1[is,:,:] .= ave.modwdh * aveH.hu1[is,:,:]
        aveWH.modwdhu3[is,:,:] .= ave.modwdh * aveH.hu3[is,:,:]
        aveWH.modwdhu3ht1[is,:,:] .= ave.modwdh * aveH.hu3ht1[is,:,:]
        aveWH.modwdhu3ht3[is,:,:] .= ave.modwdh * aveH.hu3ht3[is,:,:]
        aveWH.modwdhu33[is,:,:]   .= ave.modwdh * aveH.hu33[is,:,:]
        aveWH.modwdhu33ht1[is,:,:].= ave.modwdh * aveH.hu33ht1[is,:,:]
        aveWH.modwdhu33ht3[is,:,:].= ave.modwdh * aveH.hu33ht3[is,:,:]
        
        if(use_bper_in)
            aveWH.wdhp1b0[is,:,:]  .= ave.wdh * aveH.hp1b0[is,:,:]
            aveWH.wdhr11b0[is,:,:] .= ave.wdh * aveH.hr11b0[is,:,:]
            aveWH.wdhr13b0[is,:,:] .= ave.wdh * aveH.hr13b0[is,:,:]
            if(vpar_model_in==0)
                aveWH.wdhp1bp[is,:,:]  .= ave.wdh * aveH.hp1bp[is,:,:]
                aveWH.wdhr11bp[is,:,:] .= ave.wdh * aveH.hr11bp[is,:,:]
                aveWH.wdhr13bp[is,:,:] .= ave.wdh * aveH.hr13bp[is,:,:]
            end
        end
    end
end





#***************************************************************
#   compute the products ave_kpar*ave_h and ave_modkpar*ave_h
#***************************************************************
function kpar_h!(inputs::InputTJLF{T},ave::Ave{T},aveH::AveH{T},aveKH::AveKH) where T<:Real

    use_bper_in = inputs.USE_BPER
    vpar_model_in = inputs.VPAR_MODEL
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    for is = ns0:ns
        aveKH.kparhnp0[is,:,:]  .= aveH.hnp0[is,:,:] * ave.kpar
        aveKH.kparhp1p0[is,:,:] .= aveH.hp1p0[is,:,:] * ave.kpar
        aveKH.kparhp3p0[is,:,:] .= aveH.hp3p0[is,:,:] * ave.kpar
        aveKH.kparhu1[is,:,:]   .= aveH.hu1[is,:,:] * ave.kpar_eff[is,:,:]
        aveKH.kparhu3[is,:,:]   .= aveH.hu3[is,:,:] * ave.kpar_eff[is,:,:]
        aveKH.kparht1[is,:,:]   .= ave.kpar_eff[is,:,:] * aveH.ht1[is,:,:]
        aveKH.kparht3[is,:,:]   .= ave.kpar_eff[is,:,:] * aveH.ht3[is,:,:]
        aveKH.modkparhu1[is,:,:].= ave.modkpar_eff[is,:,:] * aveH.hu1[is,:,:]
        aveKH.modkparhu3[is,:,:].= ave.modkpar_eff[is,:,:] * aveH.hu3[is,:,:]
        if(use_bper_in)
            aveKH.kparhp1b0[is,:,:]  .= aveH.hp1b0[is,:,:] * ave.kpar
            aveKH.kparhr11b0[is,:,:] .= aveH.hr11b0[is,:,:] * ave.kpar
            aveKH.kparhr13b0[is,:,:] .= aveH.hr13b0[is,:,:] * ave.kpar
            if(vpar_model_in==0)
                aveKH.kparhnbp[is,:,:]  .= aveH.hnbp[is,:,:] * ave.kpar
                aveKH.kparhp3bp[is,:,:] .= aveH.hp3bp[is,:,:] * ave.kpar
                aveKH.kparhp1bp[is,:,:] .= aveH.hp1bp[is,:,:] * ave.kpar
                aveKH.kparhr11bp[is,:,:].= aveH.hr11bp[is,:,:] * ave.kpar
                aveKH.kparhr13bp[is,:,:].= aveH.hr13bp[is,:,:] * ave.kpar
            end
        end
    end
end


#***************************************************************
#   compute the g_ratios
#***************************************************************
function g_ratios!(inputs::InputTJLF{T}, aveG::AveG{T}) where T<:Real
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    # compute matrix ratios
    for is = ns0:ns
        gninv = inv(aveG.gn[is,:,:])
        gp1inv = inv(aveG.gp1[is,:,:])
        gp3inv = inv(aveG.gp3[is,:,:])
    
        aveG.gt1[is,:,:]  .= aveG.gp1[is,:,:] *gninv
        aveG.gt3[is,:,:]  .= aveG.gp3[is,:,:] *gninv
        aveG.gu1[is,:,:]  .= aveG.gr11[is,:,:]*gp1inv
        aveG.gu3[is,:,:]  .= aveG.gr13[is,:,:]*gp1inv
        aveG.gu33[is,:,:] .= aveG.gr33[is,:,:]*gp3inv
        aveG.gu3gt1[is,:,:]  .= aveG.gu3[is,:,:]*aveG.gt1[is,:,:]
        aveG.gu3gt3[is,:,:]  .= aveG.gu3[is,:,:]*aveG.gt3[is,:,:]
        aveG.gu33gt1[is,:,:] .= aveG.gu33[is,:,:]*aveG.gt1[is,:,:]
        aveG.gu33gt3[is,:,:] .= aveG.gu33[is,:,:]*aveG.gt3[is,:,:]
    end

end

#***************************************************************
#   compute the average g- bessel functions
#***************************************************************
function ave_gp0!(inputs::InputTJLF{T},ave::Ave{T},aveG::AveG{T}) where T<:Real
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    for is = ns0:ns
        aveG.gnp0[is,:,:] .= aveG.gn[is,:,:]*ave.p0inv
        aveG.gp1p0[is,:,:] .= aveG.gp1[is,:,:]*ave.p0inv
        aveG.gp3p0[is,:,:] .= aveG.gp3[is,:,:]*ave.p0inv
        aveG.gr11p0[is,:,:] .= aveG.gr11[is,:,:]*ave.p0inv
        aveG.gr13p0[is,:,:] .= aveG.gr13[is,:,:]*ave.p0inv
        aveG.gr33p0[is,:,:] .= aveG.gr33[is,:,:]*ave.p0inv

        aveG.c_tor_par_gp1p0[is,:,:] .= ave.c_tor_par*aveG.gp1p0[is,:,:]
        aveG.c_tor_par_gr11p0[is,:,:] .= ave.c_tor_par*aveG.gr11p0[is,:,:]        
        aveG.c_tor_par_gr13p0[is,:,:] .= ave.c_tor_par*aveG.gr13p0[is,:,:]      
    end
end

#***************************************************************
#   compute the products ave_g*ave_b0inv
#***************************************************************
function ave_gb0!(inputs::InputTJLF{T},ave::Ave{T},aveG::AveG{T}) where T<:Real
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    for is = ns0:ns
        aveG.gnb0[is,:,:]    = aveG.gn[is,:,:]  *ave.b0inv
        aveG.gp1b0[is,:,:]   = aveG.gp1[is,:,:] *ave.b0inv
        aveG.gp3b0[is,:,:]   = aveG.gp3[is,:,:] *ave.b0inv
        aveG.gr11b0[is,:,:]  = aveG.gr11[is,:,:]*ave.b0inv
        aveG.gr13b0[is,:,:]  = aveG.gr13[is,:,:]*ave.b0inv
        aveG.gr33b0[is,:,:]  = aveG.gr33[is,:,:]*ave.b0inv
        aveG.gw113b0[is,:,:] = aveG.gw113[is,:,:]*ave.b0inv
        aveG.gw133b0[is,:,:] = aveG.gw133[is,:,:]*ave.b0inv
        aveG.gw333b0[is,:,:] = aveG.gw333[is,:,:]*ave.b0inv
    end
end

#***************************************************************
#   compute the products ave_g*ave_bpinv
#***************************************************************
function ave_gbp!(inputs::InputTJLF{T},ave::Ave{T},aveG::AveG{T}) where T<:Real
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    for is = ns0:ns
        aveG.gnbp[is,:,:]    = aveG.gn[is,:,:]*ave.bpinv
        aveG.gp1bp[is,:,:]   = aveG.gp1[is,:,:]*ave.bpinv
        aveG.gp3bp[is,:,:]   = aveG.gp3[is,:,:]*ave.bpinv
        aveG.gr11bp[is,:,:]  = aveG.gr11[is,:,:]*ave.bpinv
        aveG.gr13bp[is,:,:]  = aveG.gr13[is,:,:]*ave.bpinv
        aveG.gr33bp[is,:,:]  = aveG.gr33[is,:,:]*ave.bpinv
        aveG.gw113bp[is,:,:] = aveG.gw113[is,:,:]*ave.bpinv
        aveG.gw133bp[is,:,:] = aveG.gw133[is,:,:]*ave.bpinv
        aveG.gw333bp[is,:,:] = aveG.gw333[is,:,:]*ave.bpinv
    end
end

#***************************************************************
#   compute the products ave_wd*ave_g
#***************************************************************
function wd_g!(inputs::InputTJLF{T},ave::Ave{T},aveG::AveG{T},aveWG::AveWG{T}) where T<:Real

    use_bper_in = inputs.USE_BPER
    vpar_model_in = inputs.VPAR_MODEL
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    for is = ns0:ns
        aveWG.wdgp1p0[is,:,:] = ave.wdg*aveG.gp1p0[is,:,:]
        aveWG.wdgr11p0[is,:,:] = ave.wdg*aveG.gr11p0[is,:,:]
        aveWG.wdgr13p0[is,:,:] = ave.wdg*aveG.gr13p0[is,:,:]
        aveWG.wdgu1[is,:,:] = ave.wdg*aveG.gu1[is,:,:]
        aveWG.wdgu3[is,:,:] = ave.wdg*aveG.gu3[is,:,:]
        aveWG.wdgu33[is,:,:] = ave.wdg*aveG.gu33[is,:,:]
        aveWG.wdgt1[is,:,:]= ave.wdg*aveG.gt1[is,:,:]
        aveWG.wdgt3[is,:,:]= ave.wdg*aveG.gt3[is,:,:]
        aveWG.wdgu3gt1[is,:,:]= ave.wdg*aveG.gu3gt1[is,:,:]
        aveWG.wdgu3gt3[is,:,:]= ave.wdg*aveG.gu3gt3[is,:,:]
        aveWG.wdgu33gt1[is,:,:]= ave.wdg*aveG.gu33gt1[is,:,:]
        aveWG.wdgu33gt3[is,:,:]= ave.wdg*aveG.gu33gt3[is,:,:]
        aveWG.modwdgu1[is,:,:]= ave.modwdg*aveG.gu1[is,:,:]
        aveWG.modwdgu3[is,:,:]= ave.modwdg*aveG.gu3[is,:,:]
        aveWG.modwdgu33[is,:,:]= ave.modwdg*aveG.gu33[is,:,:]
        aveWG.modwdgt1[is,:,:]= ave.modwdg*aveG.gt1[is,:,:]
        aveWG.modwdgt3[is,:,:]= ave.modwdg*aveG.gt3[is,:,:]
        aveWG.modwdgu3gt1[is,:,:]= ave.modwdg*aveG.gu3gt1[is,:,:]
        aveWG.modwdgu3gt3[is,:,:]= ave.modwdg*aveG.gu3gt3[is,:,:]
        aveWG.modwdgu33gt1[is,:,:]= ave.modwdg*aveG.gu33gt1[is,:,:]
        aveWG.modwdgu33gt3[is,:,:]= ave.modwdg*aveG.gu33gt3[is,:,:]
        if(use_bper_in)
            aveWG.wdgp1b0[is,:,:]   =     ave.wdg*aveG.gp1b0[is,:,:]
            aveWG.wdgr11b0[is,:,:] = ave.wdg*aveG.gr11b0[is,:,:]
            aveWG.wdgr13b0[is,:,:] =  ave.wdg*aveG.gr13b0[is,:,:]
            if(vpar_model_in==0)
                aveWG.wdgp1bp[is,:,:] = ave.wdg*aveG.gp1bp[is,:,:]
                aveWG.wdgr11bp[is,:,:] = ave.wdg*aveG.gr11bp[is,:,:]
                aveWG.wdgr13bp[is,:,:] = ave.wdg*aveG.gr13bp[is,:,:]
            end
        end
    end
end


#***************************************************************
#   compute the products ave_kpar*ave_h and ave_modkpar*ave_h
#***************************************************************
function kpar_g!(inputs::InputTJLF{T},ave::Ave{T},aveG::AveG{T},aveKG::AveKG) where T<:Real

    use_bper_in = inputs.USE_BPER
    vpar_model_in = inputs.VPAR_MODEL
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    for is = ns0:ns
        aveKG.kpargnp0[is,:,:]  .= aveG.gnp0[is,:,:] * ave.kpar
        aveKG.kpargp1p0[is,:,:] .= aveG.gp1p0[is,:,:] * ave.kpar
        aveKG.kpargp3p0[is,:,:] .= aveG.gp3p0[is,:,:] * ave.kpar
        aveKG.kpargu1[is,:,:]   .= aveG.gu1[is,:,:] * ave.kpar_eff[is,:,:]
        aveKG.kpargu3[is,:,:]   .= aveG.gu3[is,:,:] * ave.kpar_eff[is,:,:]
        aveKG.kpargt1[is,:,:]   .= ave.kpar_eff[is,:,:] * aveG.gt1[is,:,:]
        aveKG.kpargt3[is,:,:]   .= ave.kpar_eff[is,:,:] * aveG.gt3[is,:,:]
        aveKG.modkpargu1[is,:,:].= ave.modkpar_eff[is,:,:] * aveG.gu1[is,:,:]
        aveKG.modkpargu3[is,:,:].= ave.modkpar_eff[is,:,:] * aveG.gu3[is,:,:]
        if(use_bper_in)
            aveKG.kpargp1b0[is,:,:]  .= aveG.gp1b0[is,:,:] * ave.kpar
            aveKG.kpargr11b0[is,:,:] .= aveG.gr11b0[is,:,:] * ave.kpar
            aveKG.kpargr13b0[is,:,:] .= aveG.gr13b0[is,:,:] * ave.kpar
            if(vpar_model_in==0)
                aveKG.kpargnbp[is,:,:]  .= aveG.gnbp[is,:,:] * ave.kpar
                aveKG.kpargp3bp[is,:,:] .= aveG.gp3bp[is,:,:] * ave.kpar
                aveKG.kpargp1bp[is,:,:] .= aveG.gp1bp[is,:,:] * ave.kpar
                aveKG.kpargr11bp[is,:,:].= aveG.gr11bp[is,:,:] * ave.kpar
                aveKG.kpargr13bp[is,:,:].= aveG.gr13bp[is,:,:] * ave.kpar
            end
        end
    end
end
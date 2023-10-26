function get_matrix(inputs::InputTJLF{T}, outputGeo::OutputGeometry{T}, outputHermite::OutputHermite{T},
                    ky::T,
                    nbasis::Int) where T<:Real

    ns::Int = inputs.NS

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


    hxn = zeros(Float64, ns, nx)
    hxp1 = zeros(Float64, ns, nx)
    hxp3 = zeros(Float64, ns, nx)
    hxr11 = zeros(Float64, ns, nx)
    hxr13 = zeros(Float64, ns, nx)
    hxr33 = zeros(Float64, ns, nx)
    hxw113 = zeros(Float64, ns, nx)
    hxw133 = zeros(Float64, ns, nx)
    hxw333 = zeros(Float64, ns, nx)
    nroot = 15 ###### hardcoded
    if(nroot>6)
        gxn = zeros(Float64, ns, nx)
        gxp1 = zeros(Float64, ns, nx)
        gxp3 = zeros(Float64, ns, nx)
        gxr11 = zeros(Float64, ns, nx)
        gxr13 = zeros(Float64, ns, nx)
        gxr33 = zeros(Float64, ns, nx)
        gxw113 = zeros(Float64, ns, nx)
        gxw133 = zeros(Float64, ns, nx)
        gxw333 = zeros(Float64, ns, nx)
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
        end
    end

    aveH.hn[abs.(aveH.hn) .< zero_cut] .= 0
    aveH.hp1[abs.(aveH.hp1) .< zero_cut] .= 0
    aveH.hp3[abs.(aveH.hp3) .< zero_cut] .= 0
    aveH.hr11[abs.(aveH.hr11) .< zero_cut] .= 0
    aveH.hr13[abs.(aveH.hr13) .< zero_cut] .= 0
    aveH.hr33[abs.(aveH.hr33) .< zero_cut] .= 0
    aveH.hw113[abs.(aveH.hw113) .< zero_cut] .= 0
    aveH.hw133[abs.(aveH.hw133) .< zero_cut] .= 0
    aveH.hw333[abs.(aveH.hw333) .< zero_cut] .= 0

    if(nroot>6)
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

    ave.wdh[abs.(ave.wdh) .< zero_cut] .= 0.0
    ave.wdg[abs.(ave.wdg) .< zero_cut] .= 0.0
    ave.b0[abs.(ave.b0) .< zero_cut] .= 0.0
    ave.lnB[abs.(ave.lnB) .< zero_cut] .= 0.0
    ave.p0inv[abs.(ave.p0inv) .< zero_cut] .= 0.0
    ave.p0[abs.(ave.p0) .< zero_cut] .= 0.0
    ave.kx[abs.(ave.kx) .< zero_cut] .= 0.0
    ave.c_tor_par[abs.(ave.c_tor_par) .< zero_cut] .= 0.0
    ave.c_tor_per[abs.(ave.c_tor_per) .< zero_cut] .= 0.0
    ave.c_par_par[abs.(ave.c_par_par) .< zero_cut] .= 0.0

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


function mult1!(C, A, B, Ctmp, Atmp, is)
    Atmp .= @view A[is,:,:]
    mul!(Ctmp, Atmp, B)
    C[is,:,:] .= Ctmp
end

function mult2!(C, A, B, Ctmp, Btmp, is)
    Btmp .= @view B[is,:,:]
    mul!(Ctmp, A, Btmp)
    C[is,:,:] .= Ctmp
end

function mult3!(C, A, B, Ctmp, Atmp, Btmp, is)
    Atmp .= @view A[is,:,:]
    Btmp .= @view B[is,:,:]
    mul!(Ctmp, Atmp, Btmp)
    C[is,:,:] .= Ctmp
end


#***************************************************************
#   compute the h_ratios
#***************************************************************
function h_ratios!(inputs::InputTJLF{T}, aveH::AveH{T}) where T<:Real

    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    _, nbasis, _ = size(aveH.hn)
    Ctmp = zeros(eltype(aveH.ht1), nbasis, nbasis)
    Atmp = zeros(eltype(aveH.hp1), nbasis, nbasis)
    Btmp = zeros(eltype(aveH.hn), nbasis, nbasis)

    # compute matrix ratios
    for is = ns0:ns
        hninv = inv(aveH.hn[is,:,:])
        hp1inv = inv(aveH.hp1[is,:,:])
        hp3inv = inv(aveH.hp3[is,:,:])

        mult1!(aveH.ht1, aveH.hp1, hninv, Ctmp, Atmp, is)
        mult1!(aveH.ht3, aveH.hp3, hninv, Ctmp, Atmp, is)
        mult1!(aveH.hu1, aveH.hr11, hp1inv, Ctmp, Atmp, is)
        mult1!(aveH.hu3, aveH.hr13, hp1inv, Ctmp, Atmp, is)
        mult1!(aveH.hu33, aveH.hr33, hp3inv, Ctmp, Atmp, is)

        mult3!(aveH.hu3ht1, aveH.hu3, aveH.ht1, Ctmp, Atmp, Btmp, is)
        mult3!(aveH.hu33ht1, aveH.hu33, aveH.ht1, Ctmp, Atmp, Btmp, is)
        mult3!(aveH.hu3ht3, aveH.hu3, aveH.ht3, Ctmp, Atmp, Btmp, is)
        mult3!(aveH.hu33ht3, aveH.hu33, aveH.ht3, Ctmp, Atmp, Btmp, is)
    end

end



#***************************************************************
#   compute the products ave_h*ave_p0inv
#***************************************************************
function ave_hp0!(inputs::InputTJLF{T},ave::Ave{T},aveH::AveH{T}) where T<:Real

    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    _, nbasis, _ = size(aveH.hn)
    Ctmp = zeros(eltype(aveH.hnp0), nbasis, nbasis)
    Atmp = zeros(eltype(aveH.hn), nbasis, nbasis)
    Btmp = zeros(eltype(aveH.hp1p0), nbasis, nbasis)

    # compute matrix products
    for is = ns0:ns
        mult1!(aveH.hnp0, aveH.hn, ave.p0inv, Ctmp, Atmp, is)
        mult1!(aveH.hp1p0, aveH.hp1, ave.p0inv, Ctmp, Atmp, is)
        mult1!(aveH.hp3p0, aveH.hp3, ave.p0inv, Ctmp, Atmp, is)
        mult1!(aveH.hr11p0, aveH.hr11, ave.p0inv, Ctmp, Atmp, is)
        mult1!(aveH.hr13p0, aveH.hr13, ave.p0inv, Ctmp, Atmp, is)
        mult1!(aveH.hr33p0, aveH.hr33, ave.p0inv, Ctmp, Atmp, is)

        mult2!(aveH.c_tor_par_hp1p0, ave.c_tor_par, aveH.hp1p0, Ctmp, Btmp, is)
        mult2!(aveH.c_tor_par_hr11p0, ave.c_tor_par, aveH.hr11p0, Ctmp, Btmp, is)
        mult2!(aveH.c_tor_par_hr13p0, ave.c_tor_par, aveH.hr13p0, Ctmp, Btmp, is)
    end

end


#***************************************************************
#   compute the products ave_h*ave_b0inv
#***************************************************************
function ave_hb0!(inputs::InputTJLF{T},ave::Ave{T},aveH::AveH{T}) where T<:Real

    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    _, nbasis, _ = size(aveH.hn)
    Ctmp = zeros(eltype(aveH.hnb0), nbasis, nbasis)
    Atmp = zeros(eltype(aveH.hn), nbasis, nbasis)

    # compute matrix products
    for is = ns0:ns
        mult1!(aveH.hnb0, aveH.hn, ave.b0inv, Ctmp, Atmp, is)
        mult1!(aveH.hp1b0, aveH.hp1, ave.b0inv, Ctmp, Atmp, is)
        mult1!(aveH.hp3b0, aveH.hp3, ave.b0inv, Ctmp, Atmp, is)
        mult1!(aveH.hr11b0, aveH.hr11, ave.b0inv, Ctmp, Atmp, is)
        mult1!(aveH.hr13b0, aveH.hr13, ave.b0inv, Ctmp, Atmp, is)
        mult1!(aveH.hr33b0, aveH.hr33, ave.b0inv, Ctmp, Atmp, is)
        mult1!(aveH.hw113b0, aveH.hw113, ave.b0inv, Ctmp, Atmp, is)
        mult1!(aveH.hw133b0, aveH.hw133, ave.b0inv, Ctmp, Atmp, is)
        mult1!(aveH.hw333b0, aveH.hw333, ave.b0inv, Ctmp, Atmp, is)
    end
end



#***************************************************************
#   compute the products ave_h*ave_bpinv
#***************************************************************

function ave_hbp!(inputs::InputTJLF{T},ave::Ave{T},aveH::AveH{T}) where T<:Real

    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    _, nbasis, _ = size(aveH.hnbp)
    Ctmp = zeros(eltype(aveH.hnbp), nbasis, nbasis)
    Atmp = zeros(eltype(aveH.hn), nbasis, nbasis)
    for is = ns0:ns
        mult1!(aveH.hnbp, aveH.hn, ave.bpinv, Ctmp, Atmp, is)
        mult1!(aveH.hp1bp, aveH.hp1, ave.bpinv, Ctmp, Atmp, is)
        mult1!(aveH.hp3bp, aveH.hp3, ave.bpinv, Ctmp, Atmp, is)
        mult1!(aveH.hr11bp, aveH.hr11, ave.bpinv, Ctmp, Atmp, is)
        mult1!(aveH.hr13bp, aveH.hr13, ave.bpinv, Ctmp, Atmp, is)
        mult1!(aveH.hr33bp, aveH.hr33, ave.bpinv, Ctmp, Atmp, is)
        mult1!(aveH.hw113bp, aveH.hw113, ave.bpinv, Ctmp, Atmp, is)
        mult1!(aveH.hw133bp, aveH.hw133, ave.bpinv, Ctmp, Atmp, is)
        mult1!(aveH.hw333bp, aveH.hw333, ave.bpinv, Ctmp, Atmp, is)
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

    _, nbasis, _ = size(aveWH.wdhp1p0)
    Ctmp = zeros(eltype(aveWH.wdhp1p0), nbasis, nbasis)
    Btmp = zeros(eltype(aveH.hr11p0), nbasis, nbasis)

    for is = ns0:ns

        mult2!(aveWH.wdhp1p0, ave.wdh, aveH.hr11p0, Ctmp, Btmp, is)
        mult2!(aveWH.wdhr11p0, ave.wdh, aveH.hr11p0, Ctmp, Btmp, is)
        mult2!(aveWH.wdhr13p0, ave.wdh, aveH.hr13p0, Ctmp, Btmp, is)
        mult2!(aveWH.wdht1, ave.wdh, aveH.ht1, Ctmp, Btmp, is)
        mult2!(aveWH.wdht3, ave.wdh, aveH.ht3, Ctmp, Btmp, is)
        mult2!(aveWH.wdhu1, ave.wdh, aveH.hu1, Ctmp, Btmp, is)
        mult2!(aveWH.wdhu3, ave.wdh, aveH.hu3, Ctmp, Btmp, is)
        mult2!(aveWH.wdhu3ht1, ave.wdh, aveH.hu3ht1, Ctmp, Btmp, is)
        mult2!(aveWH.wdhu3ht3, ave.wdh, aveH.hu3ht3, Ctmp, Btmp, is)
        mult2!(aveWH.wdhu33, ave.wdh, aveH.hu33, Ctmp, Btmp, is)
        mult2!(aveWH.wdhu33ht1, ave.wdh, aveH.hu33ht1, Ctmp, Btmp, is)
        mult2!(aveWH.wdhu33ht3, ave.wdh, aveH.hu33ht3, Ctmp, Btmp, is)
        mult2!(aveWH.modwdht1, ave.modwdh, aveH.ht1, Ctmp, Btmp, is)
        mult2!(aveWH.modwdht3, ave.modwdh, aveH.ht3, Ctmp, Btmp, is)
        mult2!(aveWH.modwdhu1, ave.modwdh, aveH.hu1, Ctmp, Btmp, is)
        mult2!(aveWH.modwdhu3, ave.modwdh, aveH.hu3, Ctmp, Btmp, is)
        mult2!(aveWH.modwdhu3ht1, ave.modwdh, aveH.hu3ht1, Ctmp, Btmp, is)
        mult2!(aveWH.modwdhu3ht3, ave.modwdh, aveH.hu3ht3, Ctmp, Btmp, is)
        mult2!(aveWH.modwdhu33, ave.modwdh, aveH.hu33, Ctmp, Btmp, is)
        mult2!(aveWH.modwdhu33ht1, ave.modwdh, aveH.hu33ht1, Ctmp, Btmp, is)
        mult2!(aveWH.modwdhu33ht3, ave.modwdh, aveH.hu33ht3, Ctmp, Btmp, is)

        if(use_bper_in)
            mult2!(aveWH.wdhp1b0, ave.wdh, aveH.hp1b0, Ctmp, Btmp, is)
            mult2!(aveWH.wdhr11b0, ave.wdh, aveH.hr11b0, Ctmp, Btmp, is)
            mult2!(aveWH.wdhr13b0, ave.wdh, aveH.hr13b0, Ctmp, Btmp, is)
            if(vpar_model_in==0)
                mult2!(aveWH.wdhp1bp, ave.wdh, aveH.hp1bp, Ctmp, Btmp, is)
                mult2!(aveWH.wdhr11bp, ave.wdh, aveH.hr11bp, Ctmp, Btmp, is)
                mult2!(aveWH.wdhr13bp, ave.wdh, aveH.hr13bp, Ctmp, Btmp, is)
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

    _, nbasis, _ = size(aveH.hnp0)
    Ctmp = zeros(eltype(aveKH.kparhnp0), nbasis, nbasis)
    Atmp = zeros(eltype(aveH.hnp0), nbasis, nbasis)
    Btmp = zeros(eltype(ave.kpar_eff), nbasis, nbasis)

    for is = ns0:ns
        mult1!(aveKH.kparhnp0, aveH.hnp0, ave.kpar, Ctmp, Atmp, is)
        mult1!(aveKH.kparhp1p0, aveH.hp1p0, ave.kpar, Ctmp, Atmp, is)
        mult1!(aveKH.kparhp3p0, aveH.hp3p0, ave.kpar, Ctmp, Atmp, is)
        mult3!(aveKH.kparhu1, aveH.hu1, ave.kpar_eff, Ctmp, Atmp, Btmp, is)
        mult3!(aveKH.kparhu3, aveH.hu3, ave.kpar_eff, Ctmp, Atmp, Btmp, is)
        mult3!(aveKH.kparht1, ave.kpar_eff, aveH.ht1, Ctmp, Btmp, Atmp, is)
        mult3!(aveKH.kparht3, ave.kpar_eff, aveH.ht3, Ctmp, Btmp, Atmp, is)
        mult3!(aveKH.modkparhu1, ave.modkpar_eff, aveH.hu1, Ctmp, Btmp, Atmp, is)
        mult3!(aveKH.modkparhu3, ave.modkpar_eff, aveH.hu3, Ctmp, Btmp, Atmp, is)
        if(use_bper_in)
            mult1!(aveKH.kparhp1b0, aveH.hp1b0, ave.kpar, Ctmp, Atmp, is)
            mult1!(aveKH.kparhr11b0, aveH.hr11b0, ave.kpar, Ctmp, Atmp, is)
            mult1!(aveKH.kparhr13b0, aveH.hr13b0, ave.kpar, Ctmp, Atmp, is)
            if(vpar_model_in==0)
                mult1!(aveKH.kparhnbp, aveH.hnbp, ave.kpar, Ctmp, Atmp, is)
                mult1!(aveKH.kparhp3bp, aveH.hp3bp, ave.kpar, Ctmp, Atmp, is)
                mult1!(aveKH.kparhp1bp, aveH.hp1bp, ave.kpar, Ctmp, Atmp, is)
                mult1!(aveKH.kparhr11bp, aveH.hr11bp, ave.kpar, Ctmp, Atmp, is)
                mult1!(aveKH.kparhr13bp, aveH.hr13bp, ave.kpar, Ctmp, Atmp, is)
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

    _, nbasis, _ = size(aveG.gn)
    Ctmp = zeros(eltype(aveG.gt1), nbasis, nbasis)
    Atmp = zeros(eltype(aveG.gp1), nbasis, nbasis)
    Btmp = zeros(eltype(aveG.gt1), nbasis, nbasis)

    # compute matrix ratios
    for is = ns0:ns
        gninv = inv(aveG.gn[is,:,:])
        gp1inv = inv(aveG.gp1[is,:,:])
        gp3inv = inv(aveG.gp3[is,:,:])

        mult1!(aveG.gt1, aveG.gp1, gninv, Ctmp, Atmp, is)
        mult1!(aveG.gt3, aveG.gp3, gninv, Ctmp, Atmp, is)
        mult1!(aveG.gu1, aveG.gr11, gp1inv, Ctmp, Atmp, is)
        mult1!(aveG.gu3, aveG.gr13, gp1inv, Ctmp, Atmp, is)
        mult1!(aveG.gu33, aveG.gr33, gp3inv, Ctmp, Atmp, is)

        mult3!(aveG.gu3gt1, aveG.gu3, aveG.gt1, Ctmp, Atmp, Btmp, is)
        mult3!(aveG.gu3gt3, aveG.gu3, aveG.gt3, Ctmp, Atmp, Btmp, is)
        mult3!(aveG.gu33gt1, aveG.gu33, aveG.gt1, Ctmp, Atmp, Btmp, is)
        mult3!(aveG.gu33gt3, aveG.gu33, aveG.gt3, Ctmp, Atmp, Btmp, is)
    end

end

#***************************************************************
#   compute the average g- bessel functions
#***************************************************************
function ave_gp0!(inputs::InputTJLF{T},ave::Ave{T},aveG::AveG{T}) where T<:Real
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    _, nbasis, _ = size(aveG.gnp0)
    Ctmp = zeros(eltype(aveG.gnp0), nbasis, nbasis)
    Atmp = zeros(eltype(aveG.gn), nbasis, nbasis)
    Btmp = zeros(eltype(aveG.gp1p0), nbasis, nbasis)

    for is = ns0:ns
        mult1!(aveG.gnp0, aveG.gn, ave.p0inv, Ctmp, Atmp, is)
        mult1!(aveG.gp1p0, aveG.gp1, ave.p0inv, Ctmp, Atmp, is)
        mult1!(aveG.gp3p0, aveG.gp3, ave.p0inv, Ctmp, Atmp, is)
        mult1!(aveG.gr11p0, aveG.gr11, ave.p0inv, Ctmp, Atmp, is)
        mult1!(aveG.gr13p0, aveG.gr13, ave.p0inv, Ctmp, Atmp, is)
        mult1!(aveG.gr33p0, aveG.gr33, ave.p0inv, Ctmp, Atmp, is)

        mult2!(aveG.c_tor_par_gp1p0, ave.c_tor_par, aveG.gp1p0, Ctmp, Btmp, is)
        mult2!(aveG.c_tor_par_gr11p0, ave.c_tor_par, aveG.gr11p0, Ctmp, Btmp, is)
        mult2!(aveG.c_tor_par_gr13p0, ave.c_tor_par, aveG.gr13p0, Ctmp, Btmp, is)
    end
end

#***************************************************************
#   compute the products ave_g*ave_b0inv
#***************************************************************
function ave_gb0!(inputs::InputTJLF{T},ave::Ave{T},aveG::AveG{T}) where T<:Real
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    _, nbasis, _ = size(aveG.gnb0)
    Ctmp = zeros(eltype(aveG.gnb0), nbasis, nbasis)
    Atmp = zeros(eltype(aveG.gn), nbasis, nbasis)

    for is = ns0:ns
        mult1!(aveG.gnb0, aveG.gn, ave.b0inv, Ctmp, Atmp, is)
        mult1!(aveG.gp1b0, aveG.gp1, ave.b0inv, Ctmp, Atmp, is)
        mult1!(aveG.gp3b0, aveG.gp3, ave.b0inv, Ctmp, Atmp, is)
        mult1!(aveG.gr11b0, aveG.gr11, ave.b0inv, Ctmp, Atmp, is)
        mult1!(aveG.gr13b0, aveG.gr13, ave.b0inv, Ctmp, Atmp, is)
        mult1!(aveG.gr33b0, aveG.gr33, ave.b0inv, Ctmp, Atmp, is)
        mult1!(aveG.gw113b0, aveG.gw113, ave.b0inv, Ctmp, Atmp, is)
        mult1!(aveG.gw133b0, aveG.gw133, ave.b0inv, Ctmp, Atmp, is)
        mult1!(aveG.gw333b0, aveG.gw333, ave.b0inv, Ctmp, Atmp, is)
    end
end

#***************************************************************
#   compute the products ave_g*ave_bpinv
#***************************************************************
function ave_gbp!(inputs::InputTJLF{T},ave::Ave{T},aveG::AveG{T}) where T<:Real
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    _, nbasis, _ = size(aveG.gnbp)
    Ctmp = zeros(eltype(aveG.gnbp), nbasis, nbasis)
    Atmp = zeros(eltype(aveG.gn), nbasis, nbasis)

    for is = ns0:ns
        mult1!(aveG.gnbp, aveG.gn, ave.bpinv, Ctmp, Atmp, is)
        mult1!(aveG.gp1bp, aveG.gp1, ave.bpinv, Ctmp, Atmp, is)
        mult1!(aveG.gp3bp, aveG.gp3, ave.bpinv, Ctmp, Atmp, is)
        mult1!(aveG.gr11bp, aveG.gr11, ave.bpinv, Ctmp, Atmp, is)
        mult1!(aveG.gr13bp, aveG.gr13, ave.bpinv, Ctmp, Atmp, is)
        mult1!(aveG.gr33bp, aveG.gr33, ave.bpinv, Ctmp, Atmp, is)
        mult1!(aveG.gw113bp, aveG.gw113, ave.bpinv, Ctmp, Atmp, is)
        mult1!(aveG.gw133bp, aveG.gw133, ave.bpinv, Ctmp, Atmp, is)
        mult1!(aveG.gw333bp, aveG.gw333, ave.bpinv, Ctmp, Atmp, is)
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

    _, nbasis, _ = size(aveWG.wdgp1p0)
    Ctmp = zeros(eltype(aveWG.wdgp1p0), nbasis, nbasis)
    Btmp = zeros(eltype(aveG.gp1p0), nbasis, nbasis)

    for is = ns0:ns
        mult2!(aveWG.wdgp1p0, ave.wdg, aveG.gp1p0, Ctmp, Btmp, is)
        mult2!(aveWG.wdgr11p0, ave.wdg, aveG.gr11p0, Ctmp, Btmp, is)
        mult2!(aveWG.wdgr13p0, ave.wdg, aveG.gr13p0, Ctmp, Btmp, is)
        mult2!(aveWG.wdgu1, ave.wdg, aveG.gu1, Ctmp, Btmp, is)
        mult2!(aveWG.wdgu3, ave.wdg, aveG.gu3, Ctmp, Btmp, is)
        mult2!(aveWG.wdgu33, ave.wdg, aveG.gu33, Ctmp, Btmp, is)
        mult2!(aveWG.wdgt1, ave.wdg, aveG.gt1, Ctmp, Btmp, is)
        mult2!(aveWG.wdgt3, ave.wdg, aveG.gt3, Ctmp, Btmp, is)
        mult2!(aveWG.wdgu3gt1, ave.wdg, aveG.gu3gt1, Ctmp, Btmp, is)
        mult2!(aveWG.wdgu3gt3, ave.wdg, aveG.gu3gt3, Ctmp, Btmp, is)
        mult2!(aveWG.wdgu33gt1, ave.wdg, aveG.gu33gt1, Ctmp, Btmp, is)
        mult2!(aveWG.wdgu33gt3, ave.wdg, aveG.gu33gt3, Ctmp, Btmp, is)
        mult2!(aveWG.modwdgu1, ave.modwdg, aveG.gu1, Ctmp, Btmp, is)
        mult2!(aveWG.modwdgu3, ave.modwdg, aveG.gu3, Ctmp, Btmp, is)
        mult2!(aveWG.modwdgu33, ave.modwdg, aveG.gu33, Ctmp, Btmp, is)
        mult2!(aveWG.modwdgt1, ave.modwdg, aveG.gt1, Ctmp, Btmp, is)
        mult2!(aveWG.modwdgt3, ave.modwdg, aveG.gt3, Ctmp, Btmp, is)
        mult2!(aveWG.modwdgu3gt1, ave.modwdg, aveG.gu3gt1, Ctmp, Btmp, is)
        mult2!(aveWG.modwdgu3gt3, ave.modwdg, aveG.gu3gt3, Ctmp, Btmp, is)
        mult2!(aveWG.modwdgu33gt1, ave.modwdg, aveG.gu33gt1, Ctmp, Btmp, is)
        mult2!(aveWG.modwdgu33gt3, ave.modwdg, aveG.gu33gt3, Ctmp, Btmp, is)
        if(use_bper_in)
            mult2!(aveWG.wdgp1b0, ave.wdg, aveG.gp1b0, Ctmp, Btmp, is)
            mult2!(aveWG.wdgr11b0, ave.wdg, aveG.gr11b0, Ctmp, Btmp, is)
            mult2!(aveWG.wdgr13b0, ave.wdg, aveG.gr13b0, Ctmp, Btmp, is)
            if(vpar_model_in==0)
                mult2!(aveWG.wdgp1bp, ave.wdg, aveG.gp1bp, Ctmp, Btmp, is)
                mult2!(aveWG.wdgr11bp, ave.wdg, aveG.gr11bp, Ctmp, Btmp, is)
                mult2!(aveWG.wdgr13bp, ave.wdg, aveG.gr13bp, Ctmp, Btmp, is)
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

    _, nbasis, _ = size(aveKG.kpargnp0)
    Ctmp = zeros(eltype(aveKG.kpargnp0), nbasis, nbasis)
    Atmp = zeros(eltype(aveG.gnp0), nbasis, nbasis)
    Btmp = zeros(eltype(ave.kpar_eff), nbasis, nbasis)

    for is = ns0:ns
        mult1!(aveKG.kpargnp0, aveG.gnp0, ave.kpar, Ctmp, Atmp, is)
        mult1!(aveKG.kpargp1p0, aveG.gp1p0, ave.kpar, Ctmp, Atmp, is)
        mult1!(aveKG.kpargp3p0, aveG.gp3p0, ave.kpar, Ctmp, Atmp, is)
        mult3!(aveKG.kpargu1, aveG.gu1, ave.kpar_eff, Ctmp, Atmp, Btmp, is)
        mult3!(aveKG.kpargu3, aveG.gu3, ave.kpar_eff, Ctmp, Atmp, Btmp, is)
        mult3!(aveKG.kpargt1, ave.kpar_eff, aveG.gt1, Ctmp, Btmp, Atmp, is)
        mult3!(aveKG.kpargt3, ave.kpar_eff, aveG.gt3, Ctmp, Btmp, Atmp, is)
        mult3!(aveKG.modkpargu1, ave.modkpar_eff, aveG.gu1, Ctmp, Btmp, Atmp, is)
        mult3!(aveKG.modkpargu3, ave.modkpar_eff, aveG.gu3, Ctmp, Btmp, Atmp, is)
        if(use_bper_in)
            mult1!(aveKG.kpargp1b0, aveG.gp1b0, ave.kpar, Ctmp, Atmp, is)
            mult1!(aveKG.kpargr11b0, aveG.gr11b0, ave.kpar, Ctmp, Atmp, is)
            mult1!(aveKG.kpargr13b0, aveG.gr13b0, ave.kpar, Ctmp, Atmp, is)
            if(vpar_model_in==0)
                mult1!(aveKG.kpargnbp, aveG.gnbp, ave.kpar, Ctmp, Atmp, is)
                mult1!(aveKG.kpargp3bp, aveG.gp3bp, ave.kpar, Ctmp, Atmp, is)
                mult1!(aveKG.kpargp1bp, aveG.gp1bp, ave.kpar, Ctmp, Atmp, is)
                mult1!(aveKG.kpargr11bp, aveG.gr11bp, ave.kpar, Ctmp, Atmp, is)
                mult1!(aveKG.kpargr13bp, aveG.gr13bp, ave.kpar, Ctmp, Atmp, is)
            end
        end
    end
end
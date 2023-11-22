"""
    function get_matrix(inputs::InputTJLF{T}, outputGeo::OutputGeometry{T}, outputHermite::OutputHermite{T},ky::T,nbasis::Int)

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    outputGeo::OutputGeometry{T}        - OutputGeometry struct constructed in tjlf_geometry.jl
    outputHermite::OutputHermite{T}     - OutputHermite struct constructed in tjlf_hermite.jl
    ky:T                                - value of ky
    nbasis::T                           - number of basis for matrix dimension

outputs:
    ave
    aveH
    aveWH
    aveKH
    aveG
    aveWG
    aveKG

description:
    calculates components used for the eigenmatrix with helper functions that calculate the finite larmor radius
    integral and hermite basis averages
"""
function get_matrix(inputs::InputTJLF{T}, outputGeo::OutputGeometry{T}, outputHermite::OutputHermite{T},
                    ky::T,
                    nbasis::Int, ky_index::Int) where T<:Real

    ns::Int = inputs.NS

    ave = Ave{T}(ns, nbasis)
    aveH = AveH{T}(ns, nbasis)
    aveWH = AveWH{T}(ns, nbasis)
    aveKH = AveKH(ns, nbasis)

    aveG = AveG{T}(ns, nbasis)
    aveWG = AveWG{T}(ns, nbasis)
    aveKG = AveKG(ns, nbasis)

    aveGrad = AveGrad{Float64}(ns, nbasis)
    aveGradB = AveGradB{Float64}(ns, nbasis)

    FLR_xgrid!(inputs, outputGeo, outputHermite, aveH, aveG, ky, nbasis, ky_index)
    get_ave!(inputs, outputGeo, outputHermite, ave, nbasis, ky, ky_index)

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

    if(inputs.LINSKER_FACTOR!=0.0) 
        @warn "NOT TESTED matrix.jl ln 68"
        grad_ave_h(inputs,ave,aveH,aveGrad) 
    end

    ave_hp0!(inputs, ave, aveH)
    if(inputs.BETAE>0.0)
        ave_hb0!(inputs, ave, aveH)
        if(inputs.VPAR_MODEL==0) ave_hbp!(inputs, ave, aveH) end
    end

    wd_h!(inputs, ave, aveH, aveWH)
    kpar_h!(inputs, ave, aveH, aveKH)

    if(inputs.GRADB_FACTOR!=0.0) 
        @warn "NOT TESTED matrix.jl ln 82"
        gradB_h(inputs,ave,aveH,aveGradB) 
    end


    nroot = 15 ###### hardcoded
    if(nroot>6)
        g_ratios!(inputs, aveG)
        if(inputs.LINSKER_FACTOR!=0) 
            @warn "NOT TESTED matrix.jl ln 91"
            grad_ave_g(inputs,ave,aveG,aveGrad)
        end
        ave_gp0!(inputs, ave, aveG)
        if(inputs.BETAE>0.0)
            ave_gb0!(inputs, ave, aveG)
            if(inputs.VPAR_MODEL==0) ave_gbp!(inputs, ave, aveG) end
        end
        wd_g!(inputs, ave, aveG, aveWG)
        kpar_g!(inputs, ave, aveG, aveKG)
        if(inputs.GRADB_FACTOR!=0.0)
            @warn "NOT TESTED matrix.jl ln 102"
            gradB_g(inputs,ave,aveG,aveGradB)
        end
    end

    return ave, aveH, aveWH, aveKH, aveG, aveWG, aveKG, aveGrad, aveGradB
end

function outer!(Y, x, dvec, is)
    _, M, N = size(Y)
    for j in 1:N
        for i in 1:M
            Y[is, i, j] = sum(dvec[k] * x[i, k] * conj(x[j, k]) for k in eachindex(dvec))
        end
    end
    return Y
end


#*************************************************************
#  begin computation of the FLR integrals
#*************************************************************
#  compute the FLR integrals at the hermite nodes
function FLR_xgrid!(inputs::InputTJLF{T}, outputGeo::OutputGeometry{T}, outputHermite::OutputHermite{T},
    aveH::AveH{T}, aveG::AveG{T}, ky::T, nbasis::Int, ky_index::Int) where T<:Real

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
            b2 = b ^ 2
            b25 = b ^ 2.5
            hxn[is,i]    = FLR_Hn(fth, b; b2, b25)
            hxp1[is,i]   = hxn[is,i]
            hxp3[is,i]   = FLR_dHp3(fth, b; b2, b25) + hxn[is,i]
            hxr11[is,i]  = 3.0 * hxp1[is,i]
            hxr13[is,i]  = FLR_dHr13(fth, b; b2, b25) + (5/3)*hxp1[is,i]
            hxr33[is,i]  = FLR_dHr33(fth, b; b2, b25) + (5/3)*hxp3[is,i]
            hxw113[is,i] = FLR_dHw113(fth, b; b2, b25)+ (7/3)*hxr11[is,i]
            hxw133[is,i] = FLR_dHw133(fth, b; b2, b25)+ (7/3)*hxr13[is,i]
            hxw333[is,i] = FLR_dHw333(fth, b; b2, b25)+ (7/3)*hxr33[is,i]

            #***************************************************************
            #   compute the average g- bessel functions
            #***************************************************************
            if(nroot>6)
                fts = outputGeo.fts
                ftx::Float64 = fts[is]
                ft2 = ftx^2
                gxn[is,i]    = FLR_Hn(ftx, b; b2, b25)
                gxp1[is,i]   = FLR_dHp1(ftx, b; b2, b25)  + ft2*gxn[is,i]
                gxp3[is,i]   = FLR_dHp3(ftx, b; b2, b25)  + gxn[is,i]
                gxr11[is,i]  = FLR_dHr11(ftx, b; b2, b25) + 3*ft2*gxp1[is,i]
                gxr13[is,i]  = FLR_dHr13(ftx, b; b2, b25) + (5/3)*gxp1[is,i]
                gxr33[is,i]  = FLR_dHr33(ftx, b; b2, b25) + (5/3)*gxp3[is,i]
                gxw113[is,i] = FLR_dHw113(ftx, b; b2, b25)+ (7/3)*gxr11[is,i]
                gxw133[is,i] = FLR_dHw133(ftx, b; b2, b25)+ (7/3)*gxr13[is,i]
                gxw333[is,i] = FLR_dHw333(ftx, b; b2, b25)+ (7/3)*gxr33[is,i]

            end
        end
    end
    
    @views dvec = outputHermite._dvec[:,ky_index]
    for is = ns0:ns
        @views dvec .= hxn[is,:] .* wx
        outer!(aveH.hn, h, dvec, is)
        @views dvec .= hxp1[is,:] .* wx
        outer!(aveH.hp1, h, dvec, is)
        @views dvec .= hxp3[is,:] .* wx
        outer!(aveH.hp3, h, dvec, is)
        @views dvec .= hxr11[is,:] .* wx
        outer!(aveH.hr11, h, dvec, is)
        @views dvec .= hxr13[is,:] .* wx
        outer!(aveH.hr13, h, dvec, is)
        @views dvec .= hxr33[is,:] .* wx
        outer!(aveH.hr33, h, dvec, is)
        @views dvec .= hxw113[is,:] .* wx
        outer!(aveH.hw113, h, dvec, is)
        @views dvec .= hxw133[is,:] .* wx
        outer!(aveH.hw133, h, dvec, is)
        @views dvec .= hxw333[is,:] .* wx
        outer!(aveH.hw333, h, dvec, is)

        if(nroot>6)
            @views dvec .= gxn[is,:] .* wx
            outer!(aveG.gn, h, dvec, is)
            @views dvec .= gxp1[is,:] .* wx
            outer!(aveG.gp1, h, dvec, is)
            @views dvec .= gxp3[is,:] .* wx
            outer!(aveG.gp3, h, dvec, is)
            @views dvec .= gxr11[is,:] .* wx
            outer!(aveG.gr11, h, dvec, is)
            @views dvec .= gxr13[is,:] .* wx
            outer!(aveG.gr13, h, dvec, is)
            @views dvec .= gxr33[is,:] .* wx
            outer!(aveG.gr33, h, dvec, is)
            @views dvec .= gxw113[is,:] .* wx
            outer!(aveG.gw113, h, dvec, is)
            @views dvec .= gxw133[is,:] .* wx
            outer!(aveG.gw133, h, dvec, is)
            @views dvec .= gxw333[is,:] .* wx
            outer!(aveG.gw333, h, dvec, is)
        end
    end

    cut = abs.(aveH.hn) .< zero_cut
    aveH.hn[cut] .= 0
    cut .= abs.(aveH.hp1) .< zero_cut
    aveH.hp1[cut] .= 0
    cut .= abs.(aveH.hp3) .< zero_cut
    aveH.hp3[cut] .= 0
    cut .= abs.(aveH.hr11) .< zero_cut
    aveH.hr11[cut] .= 0
    cut .= abs.(aveH.hr13) .< zero_cut
    aveH.hr13[cut] .= 0
    cut .= abs.(aveH.hr33) .< zero_cut
    aveH.hr33[cut] .= 0
    cut .= abs.(aveH.hw113) .< zero_cut
    aveH.hw113[cut] .= 0
    cut .= abs.(aveH.hw133) .< zero_cut
    aveH.hw133[cut] .= 0
    cut .= abs.(aveH.hw333) .< zero_cut
    aveH.hw333[cut] .= 0

    if(nroot>6)
        cut .= abs.(aveG.gn) .< zero_cut
        aveG.gn[cut] .= 0
        cut .= abs.(aveG.gp1) .< zero_cut
        aveG.gp1[cut] .= 0
        cut .= abs.(aveG.gp3) .< zero_cut
        aveG.gp3[cut] .= 0
        cut .= abs.(aveG.gr11) .< zero_cut
        aveG.gr11[cut] .= 0
        cut .= abs.(aveG.gr13) .< zero_cut
        aveG.gr13[cut] .= 0
        cut .= abs.(aveG.gr33) .< zero_cut
        aveG.gr33[cut] .= 0
        cut .= abs.(aveG.gw113) .< zero_cut
        aveG.gw113[cut] .= 0
        cut .= abs.(aveG.gw133) .< zero_cut
        aveG.gw133[cut] .= 0
        cut .= abs.(aveG.gw333) .< zero_cut
        aveG.gw333[cut] .= 0
    end

end




#***********************************************************
#  compute  k-independent hermite basis averages
#***********************************************************
function get_ave!(inputs::InputTJLF{T},outputGeo::OutputGeometry{T},outputHermite::OutputHermite{T},ave::Ave{T},
    nbasis::Int, ky::T, ky_index::Int) where T<:Real

    zero_cut = 1.e-12
    fts = outputGeo.fts
    ft2 = fts[1]^2

    width_in = inputs.WIDTH_SPECTRUM[ky_index]

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

    cut = abs.(ave.wdh) .< zero_cut
    ave.wdh[cut] .= 0.0
    cut .= abs.(ave.wdg) .< zero_cut
    ave.wdg[cut] .= 0.0
    cut .= abs.(ave.b0) .< zero_cut
    ave.b0[cut] .= 0.0
    cut .= abs.(ave.lnB) .< zero_cut
    ave.lnB[cut] .= 0.0
    cut .= abs.(ave.p0inv) .< zero_cut
    ave.p0inv[cut] .= 0.0
    cut .= abs.(ave.p0) .< zero_cut
    ave.p0[cut] .= 0.0
    cut .= abs.(ave.kx) .< zero_cut
    ave.kx[cut] .= 0.0
    cut .= abs.(ave.c_tor_par) .< zero_cut
    ave.c_tor_par[cut] .= 0.0
    cut .= abs.(ave.c_tor_per) .< zero_cut
    ave.c_tor_per[cut] .= 0.0
    cut .= abs.(ave.c_par_par) .< zero_cut
    ave.c_par_par[cut] .= 0.0

    ave.gradB .= (ave.kpar*ave.lnB) .- (ave.lnB*ave.kpar)

    for is = ns0:ns
        if(nbasis==1)
            ave.kpar_eff[is,1,1] = -im/√(2)
            if(vpar_model_in==1)
                @warn "I haven't tested this (matrix.jl ln 394)"
                vpar_s = inputs.ALPHA_MACH*inputs.SIGN_IT*inputs.VPAR[is]
                ave.kpar_eff[is,1,1] = ave.kpar_eff[is,1,1] + im*abs(alpha_mach_in*vpar_s)*ky*R_unit*q_unit*inputs.WIDTH_SPECTRUM[ky_index]*inputs.MASS[is]/inputs.ZS[is]
            end
        else
            for i = 1:nbasis
                for j = 1:nbasis
                    ave.kpar_eff[is,i,j] = ave.kpar[i,j]
                    if(vpar_model_in==1 && i==j)
                        @warn "I haven't tested this (matrix.jl ln 403)"
                        vpar_s = inputs.ALPHA_MACH*inputs.SIGN_IT*inputs.VPAR[is]
                        ave.kpar_eff[is,i,j] = ave.kpar_eff[is,i,j] - im*alpha_mach_in*vpar_s*ky*R_unit*q_unit*inputs.WIDTH_SPECTRUM[ky_index]*inputs.MASS[is]/inputs.ZS[is]
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
    eigen = eigen!(Symmetric(ave.wdh))
    w = eigen.values
    a = eigen.vectors
    # if abs(w)<wd_zero, make it +- wd_zero. otherwise keep the w
    w .= sign.(w) .* max.(abs.(w), wd_zero_in)
    w = Diagonal(w)
    ave.modwdh = a * abs.(w) * transpose(a)
    ave.wdh = a * w * transpose(a) 

    ### now for wdg
    eigen = eigen!(Symmetric(ave.wdg))
    w = eigen.values
    a = eigen.vectors
    # if abs(w)<wd_zero, make it +- wd_zero. otherwise keep the w
    w .= sign.(w) .* max.(abs.(w), wd_zero_in)
    w = Diagonal(w)
    ave.modwdg = a * abs.(w) * transpose(a)
    ave.wdg = a * w * transpose(a)
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
            a = Hermitian(im.*ave.kpar)
            eigen = eigen!(a)
            w = Diagonal(eigen.values)
            a = eigen.vectors
            b = a * abs.(w) * adjoint(a)
            ave.modkpar .= real.(b)
            for is = ns0:ns
                ave.modkpar_eff[is,:,:] .= b
            end

        else
            for is = ns0:ns
                a = Hermitian(im.*ave.kpar_eff[is,:,:])
                eigen = eigen!(a)
                w = Diagonal(eigen.values)
                a = eigen.vectors
                b = a * abs.(w) * adjoint(a)
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





function mult4!(C, A, B, C1tmp, C2tmp, Atmp, is)
    Atmp .= @view A[is,:,:]
    mul!(C1tmp, Atmp, B)
    mul!(C2tmp, B, Atmp)
    C[is,:,:] .= C1tmp .- C2tmp
end

function mult5!(C, A, B, C1tmp, C2tmp, Btmp, is)
    Btmp .= @view B[is,:,:]
    mul!(C1tmp, A, Btmp)
    mul!(C2tmp, Btmp, A)
    C[is,:,:] .= C1tmp .- C2tmp
end

#***************************************************************
#   compute the parallel gradient of ave_m matrix
#***************************************************************
function grad_ave_h!(inputs::InputTJLF{T},ave::Ave{T},aveH::AveH{T},aveGrad::AveGrad{T}) where T<:Real

    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    _, nbasis, _ = size(aveGrad.gradhp1)
    C1tmp = zeros(eltype(aveGrad.gradhp1), nbasis, nbasis)
    C2tmp = zeros(eltype(aveGrad.gradhp1), nbasis, nbasis)
    Atmp = zeros(eltype(aveH.hu1), nbasis, nbasis)
    Btmp = zeros(eltype(aveH.hp1), nbasis, nbasis)

    hp1inv = inv(aveH.hp1)

    for is = ns0:ns
        mult5!(aveGrad.gradhp1, ave.kpar, aveH.hp1, C1tmp, C2tmp, Btmp, is)
        mult3!(aveGrad.gradhr11, aveH.hu1, aveGrad.gradhp1, C1tmp, Atmp, Btmp, is)
        mult3!(aveGrad.gradhr13, aveH.hu3, aveGrad.gradhp1, C1tmp, Atmp, Btmp, is)

        mult3!(aveGrad.gradhp1p1, aveGrad.gradhp1, hp1inv, C1tmp, Atmp, Btmp, is)
        mult3!(aveGrad.gradhr11p1, aveGrad.gradhr11, hp1inv, C1tmp, Atmp, Btmp, is)
        mult3!(aveGrad.gradhr13p1, aveGrad.gradhr13, hp1inv, C1tmp, Atmp, Btmp, is)

        mult1!(aveGrad.gradhp1p0, aveGrad.gradhp1, ave.p0inv, C1tmp, Atmp, is)
        mult1!(aveGrad.gradhr11p0, aveGrad.gradhr11, ave.p0inv, C1tmp, Atmp, is)
        mult1!(aveGrad.gradhr13p0, aveGrad.gradhr13, ave.p0inv, C1tmp, Atmp, is)
    end 
end


#***************************************************************
#   compute the parallel gradient of ave_g matricies
#***************************************************************
function  grad_ave_g(inputs::InputTJLF{T},ave::Ave{T},aveG::AveG{T},aveGrad::AveGrad{T}) where T<:Real

    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    _, nbasis, _ = size(aveGrad.gradgp1)
    C1tmp = zeros(eltype(aveGrad.gradgp1), nbasis, nbasis)
    C2tmp = zeros(eltype(aveGrad.gradgp1), nbasis, nbasis)
    Atmp = zeros(eltype(aveG.gu1), nbasis, nbasis)
    Btmp = zeros(eltype(aveG.gp1), nbasis, nbasis)

    gp1inv = inv(aveG.gp1)

    for is = ns0:ns
        mult5!(aveGrad.gradgp1, ave.kpar, aveG.gp1, C1tmp, C2tmp, Btmp, is)
        mult3!(aveGrad.gradgr11, aveG.gu1, aveGrad.gradgp1, C1tmp, Atmp, Btmp, is)
        mult3!(aveGrad.gradgr13, aveG.gu3, aveGrad.gradgp1, C1tmp, Atmp, Btmp, is)

        mult3!(aveGrad.gradgp1p1, aveGrad.gradgp1, gp1inv, C1tmp, Atmp, Btmp, is)
        mult3!(aveGrad.gradgr11p1, aveGrad.gradgr11, gp1inv, C1tmp, Atmp, Btmp, is)
        mult3!(aveGrad.gradgr13p1, aveGrad.gradgr13, gp1inv, C1tmp, Atmp, Btmp, is)

        mult1!(aveGrad.gradgp1p0, aveGrad.gradgp1, ave.p0inv, C1tmp, Atmp, is)
        mult1!(aveGrad.gradgr11p0, aveGrad.gradgr11, ave.p0inv, C1tmp, Atmp, is)
        mult1!(aveGrad.gradgr13p0, aveGrad.gradgr13, ave.p0inv, C1tmp, Atmp, is)
    end 
end


#***************************************************************
#   compute the products gradB*h
#***************************************************************
function gradB_h(inputs::InputTJLF{T},ave::Ave{T},aveH::AveH{T},aveGradB::AveGradB{T}) where T<:Real

    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    _, nbasis, _ = size(aveGradB.gradBhp1)
    Ctmp = zeros(eltype(aveGradB.gradBhp1), nbasis, nbasis)
    Btmp = zeros(eltype(aveH.hp1p0), nbasis, nbasis)

    for is = ns0:ns
        mult2!(aveGradB.gradBhp1, ave.gradB, aveH.hp1p0, Ctmp, Btmp, is)
        mult2!(aveGradB.gradBhp3, ave.gradB, aveH.hp3p0, Ctmp, Btmp, is)
        mult2!(aveGradB.gradBhr11, ave.gradB, aveH.hr11p0, Ctmp, Btmp, is)
        mult2!(aveGradB.gradBhr13, ave.gradB, aveH.hr13p0, Ctmp, Btmp, is)
        mult2!(aveGradB.gradBhr33, ave.gradB, aveH.hr33p0, Ctmp, Btmp, is)
        mult2!(aveGradB.gradBhu1, ave.gradB, aveH.hu1, Ctmp, Btmp, is)
        mult2!(aveGradB.gradBhu3, ave.gradB, aveH.hu3, Ctmp, Btmp, is)
        mult2!(aveGradB.gradBhu33, ave.gradB, aveH.hu33, Ctmp, Btmp, is)
    end
end


#***************************************************************
#   compute the products gradB*ave_g
#***************************************************************
function gradB_g(inputs::InputTJLF{T},ave::Ave{T},aveG::AveH{T},aveGradB::AveGradB{T}) where T<:Real

    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    _, nbasis, _ = size(aveGradB.gradBgp1)
    Ctmp = zeros(eltype(aveGradB.gradBgp1), nbasis, nbasis)
    Btmp = zeros(eltype(aveG.gp1p0), nbasis, nbasis)

    for is = ns0:ns
        mult2!(aveGradB.gradBgp1, ave.gradB, aveG.gp1p0, Ctmp, Btmp, is)
        mult2!(aveGradB.gradBgp3, ave.gradB, aveG.gp3p0, Ctmp, Btmp, is)
        mult2!(aveGradB.gradBgr11, ave.gradB, aveG.gr11p0, Ctmp, Btmp, is)
        mult2!(aveGradB.gradBgr13, ave.gradB, aveG.gr13p0, Ctmp, Btmp, is)
        mult2!(aveGradB.gradBgr33, ave.gradB, aveG.gr33p0, Ctmp, Btmp, is)
        mult2!(aveGradB.gradBgu1, ave.gradB, aveG.gu1, Ctmp, Btmp, is)
        mult2!(aveGradB.gradBgu3, ave.gradB, aveG.gu3, Ctmp, Btmp, is)
        mult2!(aveGradB.gradBgu33, ave.gradB, aveG.gu33, Ctmp, Btmp, is)
    end
end
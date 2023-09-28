include("tjlf_modules.jl")


function get_matrix()

    # FLR_xgrid()
    # ave_theta()

    # if(inputs.VPAR_MODEL==0 && inputs.USE_BPER)
    #     ave_inv0(ave_bp,ave_bpinv)
    #     for i = 1:nbasis
    #         for j = 1:nbasis
    #             ave_p0inv(i,j) = 0.0
    #             ave_b0inv(i,j) = 0.0
    #             for k = 1:nbasis
    #                 ave_p0inv(i,j) = ave_p0inv(i,j)+ave_bpinv(i,k)*ave_b0(k,j)
    #                 ave_b0inv(i,j) = ave_b0inv(i,j)+ave_bpinv(i,k)*ave_p0(k,j)
    #             end
    #         end
    #     end

    # else
    #     ave_inv0(ave_p0,ave_p0inv)
    #     if(inputs.USE_BPER) ave_inv0(ave_b0,ave_b0inv) end
    # end

    # modwd()
    # modkpar()
    # ave_h()

    # ave_inv(ave_hn,ave_hninv)
    # ave_inv(ave_hp1,ave_hp1inv)
    # ave_inv(ave_hp3,ave_hp3inv)

    # h_ratios()

    # if(inputs.LINSKER_FACTOR!=0.0) grad_ave_h() end

    # ave_hp0()
    # if(inputs.BETAE >0.0)
    #     ave_hb0()
    #     if(inputs.VPAR_MODEL==0) ave_hbp() end
    # end

    # wd_h()
    # kpar_h()
       
    # if(inputs.GRADB_FACTOR!=0.0) gradB_h() end

    # if(nroot>6)
    #     ave_g()
    #     ave_inv(ave_gn,ave_gninv)
    #     ave_inv(ave_gp1,ave_gp1inv) 
    #     ave_inv(ave_gp3,ave_gp3inv)

    #     g_ratios()

    #     if(inputs.LINSKER_FACTOR!=0.0) grad_ave_g() end

    #     ave_gp0()
    #     if(inputs.BETAE>0.0)
    #         ave_gb0()
    #         if(inputs.VPAR_MODEL==0) ave_gbp() end
    #     end

    #     wd_g()
    #     kpar_g()

    #     if(inputs.GRADB_FACTOR!=0.0) gradB_g() end

    # end

    fileDirectory = "./ave_p0inv"
    lines = readlines(fileDirectory)
    nbasis = parse(Int,split(lines[1])[3])
    ave_p0inv = Vector{Float64}()
    for line in lines[2:length(lines)]
        if contains(line,"=") continue end
        push!(ave_p0inv,parse(Float64,strip(line)))
    end
    ave_p0inv = Matrix(reshape(ave_p0inv, nbasis, nbasis)')

    fileDirectory = "./ave_b0inv"
    lines = readlines(fileDirectory)
    nbasis = parse(Int,split(lines[1])[3])
    ave_b0inv = Vector{Float64}()
    for line in lines[2:length(lines)]
        if contains(line,"=") continue end
        push!(ave_b0inv,parse(Float64,strip(line)))
    end
    ave_b0inv = Matrix(reshape(ave_b0inv, nbasis, nbasis)')

    fileDirectory = "./ave_bpinv"
    lines = readlines(fileDirectory)
    nbasis = parse(Int,split(lines[1])[3])
    ave_bpinv = Vector{Float64}()
    for line in lines[2:length(lines)]
        if contains(line,"=") continue end
        push!(ave_bpinv,parse(Float64,strip(line)))
    end
    ave_bpinv = Matrix(reshape(ave_bpinv, nbasis, nbasis)')

    fileDirectory = "./ave_wdh"
    lines = readlines(fileDirectory)
    nbasis = parse(Int,split(lines[1])[3])
    ave_wdh = Vector{Float64}()
    for line in lines[2:length(lines)]
        if contains(line,"=") continue end
        push!(ave_wdh,parse(Float64,strip(line)))
    end
    ave_wdh = Matrix(reshape(ave_wdh, nbasis, nbasis)')
    
    fileDirectory = "./ave_b0"
    lines = readlines(fileDirectory)
    nbasis = parse(Int,split(lines[1])[3])
    ave_b0 = Vector{Float64}()
    for line in lines[2:length(lines)]
        if contains(line,"=") continue end
        push!(ave_b0,parse(Float64,strip(line)))
    end
    ave_b0 = Matrix(reshape(ave_b0, nbasis, nbasis)')

    fileDirectory = "./ave_kx"
    lines = readlines(fileDirectory)
    nbasis = parse(Int,split(lines[1])[3])
    ave_kx = Vector{Float64}()
    for line in lines[2:length(lines)]
        if contains(line,"=") continue end
        push!(ave_kx,parse(Float64,strip(line)))
    end
    ave_kx = Matrix(reshape(ave_kx, nbasis, nbasis)')

    fileDirectory = "./ave_c_par_par"
    lines = readlines(fileDirectory)
    nbasis = parse(Int,split(lines[1])[3])
    ave_c_par_par = Vector{Float64}()
    for line in lines[2:length(lines)]
        if contains(line,"=") continue end
        push!(ave_c_par_par,parse(Float64,strip(line)))
    end
    ave_c_par_par = Matrix(reshape(ave_c_par_par, nbasis, nbasis)')

    fileDirectory = "./ave_kpar"
    lines = readlines(fileDirectory)
    nbasis = parse(Int,split(lines[1])[3])
    ave_kpar = Vector{Float64}()
    for line in lines[2:length(lines)]
        if contains(line,"=") continue end
        push!(ave_kpar,parse(Float64,strip(line)))
    end
    ave_kpar = Matrix(reshape(ave_kpar, nbasis, nbasis)')

    fileDirectory = "./ave_c_tor_par"
    lines = readlines(fileDirectory)
    nbasis = parse(Int,split(lines[1])[3])
    ave_c_tor_par = Vector{Float64}()
    for line in lines[2:length(lines)]
        if contains(line,"=") continue end
        push!(ave_c_tor_par,parse(Float64,strip(line)))
    end
    ave_c_tor_par = Matrix(reshape(ave_c_tor_par, nbasis, nbasis)')

    fileDirectory = "./ave_c_tor_per"
    lines = readlines(fileDirectory)
    nbasis = parse(Int,split(lines[1])[3])
    ave_c_tor_per = Vector{Float64}()
    for line in lines[2:length(lines)]
        if contains(line,"=") continue end
        push!(ave_c_tor_per,parse(Float64,strip(line)))
    end
    ave_c_tor_per = Matrix(reshape(ave_c_tor_per, nbasis, nbasis)')

    fileDirectory = "./ave_hp1"
    lines = readlines(fileDirectory)
    ave_hp1 = Vector{Float64}()
    for line in lines[2:length(lines)]
        if contains(line,"=") continue end
        push!(ave_hp1,parse(Float64,strip(line)))
    end
    ave_hp1 = permutedims(reshape(ave_hp1, nbasis, nbasis, 3), (3,2,1))

    return ave_p0inv, ave_b0inv, ave_bpinv, ave_wdh, ave_b0, ave_kx, ave_c_par_par, ave_kpar, ave_c_tor_par, ave_c_tor_per, ave_hp1

end

function TJLFEP_generate_input()
    # Constants and parameters
    pi = 3.14159265359
    n_param = 0

    # File check
    filename = "input.TJLFEP.generate"
    if isfile(filename)
        open(filename, "r") do file
            sign_bt = parse(Float64, readline(file))
            sign_it = parse(Float64, readline(file))
            nr = parse(Int, readline(file))
            ns = parse(Int, readline(file))
            geometry_flag = 1
            scan_n = nr

            # Allocations
            zs = zeros(Float64, ns)
            mass = zeros(Float64, ns)

            as_min = zeros(Float64, ns)
            as_max = zeros(Float64, ns)
            taus_min = zeros(Float64, ns)
            taus_max = zeros(Float64, ns)
            rlns_min = zeros(Float64, ns)
            rlns_max = zeros(Float64, ns)
            rlts_min = zeros(Float64, ns)
            rlts_max = zeros(Float64, ns)

            # Read species data
            for i in 1:ns
                readline(file)  # Skip line
                readline(file)  # Skip line
                zs[i] = parse(Float64, readline(file))
                mass[i] = parse(Float64, readline(file))

                if i != 1
                    dum_min, dum_max = readline_values(file, n_param)
                    as_min[i] = dum_min
                    as_max[i] = dum_max
                    dum_min, dum_max = readline_values(file, n_param)
                    taus_min[i] = dum_min
                    taus_max[i] = dum_max
                else
                    as_min[i] = 1.0   # Electron density
                    as_max[i] = 1.0
                    taus_min[i] = 1.0  # Electron temperature
                    taus_max[i] = 1.0
                end

                dum_min, dum_max = readline_values(file, n_param)
                rlts_min[i] = dum_min
                rlts_max[i] = dum_max
                dum_min, dum_max = readline_values(file, n_param)
                rlns_min[i] = dum_min
                rlns_max[i] = dum_max
            end

            # Skip two lines
            readline(file)
            readline(file)

            rmin_min, rmin_max = readline_values(file, n_param)
            rmaj_min, rmaj_max = readline_values(file, n_param)
            q_min, q_max = readline_values(file, n_param)
            shear_min, shear_max = readline_values(file, n_param)
            shift_min, shift_max = readline_values(file, n_param)
            kappa_min, kappa_max = readline_values(file, n_param)
            s_kappa_min, s_kappa_max = readline_values(file, n_param)
            delta_min, delta_max = readline_values(file, n_param)
            s_delta_min, s_delta_max = readline_values(file, n_param)
            zeta_min, zeta_max = readline_values(file, n_param)
            s_zeta_min, s_zeta_max = readline_values(file, n_param)
            betae_min, betae_max = readline_values(file, n_param)

            if nr < 0
                nr = abs(nr) * n_param
                i_data_method = 1
            else
                i_data_method = 2
            end

            # Allocate arrays
            as = zeros(Float64, nr, ns)
            taus = zeros(Float64, nr, ns)
            rlns = zeros(Float64, nr, ns)
            rlts = zeros(Float64, nr, ns)

            rmin = zeros(Float64, nr)
            rmaj = zeros(Float64, nr)
            q = zeros(Float64, nr)
            shear = zeros(Float64, nr)

            q_prime = zeros(Float64, nr)
            p_prime = zeros(Float64, nr)
            shift = zeros(Float64, nr)
            kappa = zeros(Float64, nr)
            s_kappa = zeros(Float64, nr)
            delta = zeros(Float64, nr)
            s_delta = zeros(Float64, nr)
            zeta = zeros(Float64, nr)
            s_zeta = zeros(Float64, nr)

            zeff = zeros(Float64, nr)
            betae = zeros(Float64, nr)

            rho_star = zeros(Float64, nr)
            omega_TAE = zeros(Float64, nr)
            omegaGAM = zeros(Float64, nr)
            gammaE = ones(Float64, nr) * 1.0E-7
            gammap = ones(Float64, nr) * 1.0E-7

            ptot_scale = zeros(Float64, nr)
            dptotdr_scale = zeros(Float64, nr)

            i_dim = 12 + 4 * ns
            ran_arr = rand(Float64, i_dim)

            if i_data_method == 1
                println("Grid data not yet operational.")
            elseif i_data_method == 2
                for j in 1:nr
                    for i in 1:ns
                        as[j, i] = as_min[i] + ran_arr[i] * (as_max[i] - as_min[i])
                        taus[j, i] = taus_min[i] + ran_arr[i + 1] * (taus_max[i] - taus_min[i])
                        rlns[j, i] = rlns_min[i] + ran_arr[i + 2] * (rlns_max[i] - rlns_min[i])
                        rlts[j, i] = rlts_min[i] + ran_arr[i + 3] * (rlts_max[i] - rlts_min[i])
                    end

                    rmin[j] = rmin_min + ran_arr[4 * ns + 1] * (rmin_max - rmin_min)
                    rmaj[j] = rmaj_min + ran_arr[4 * ns + 2] * (rmaj_max - rmaj_min)
                    q[j] = q_min + ran_arr[4 * ns + 3] * (q_max - q_min)
                    shear[j] = shear_min + ran_arr[4 * ns + 4] * (shear_max - shear_min)
                    shift[j] = shift_min + ran_arr[4 * ns + 5] * (shift_max - shift_min)
                    kappa[j] = kappa_min + ran_arr[4 * ns + 6] * (kappa_max - kappa_min)
                    s_kappa[j] = s_kappa_min + ran_arr[4 * ns + 7] * (s_kappa_max - s_kappa_min)
                    delta[j] = delta_min + ran_arr[4 * ns + 8] * (delta_max - delta_min)
                    s_delta[j] = s_delta_min + ran_arr[4 * ns + 9] * (s_delta_max - s_delta_min)
                    zeta[j] = zeta_min + ran_arr[4 * ns + 10] * (zeta_max - zeta_min)
                    s_zeta[j] = s_zeta_min + ran_arr[4 * ns + 11] * (s_zeta_max - s_zeta_min)
                    betae[j] = betae_min + ran_arr[4 * ns + 12] * (betae_max - betae_min)
                    omega_TAE[j] = sqrt(2.0 / betae[j]) / 2.0 / q[j] / rmaj[j]
                    omegaGAM[j] = (1.0 / rmaj[j]) * sqrt(1.0 + taus[j, 2]) / (1.0 + 1.0 / (2.0 * q[j]))

                    for i in 1:ns
                        ptot_scale[j] += as[j, i] * taus[j, i]
                        dptotdr_scale[j] += (rlns[j, i] + rlts[j, i]) * as[j, i] * taus[j, i]
                        zeff[j] += as[j, i] * zs[i]
                    end
                end

                q_prime .= (q ./ rmin).^2 .* shear
                p_prime .= -1.0 * abs.(q) ./ rmin .* betae ./ (8 * pi) .* dptotdr_scale
                rho_star .= 0.001
            end
        end
    else
        println("input.TJLFEP.generate file not found")
        exit(1)
    end
end

function readline_values(file, n_count)
    line = readline(file)
    lead_char = first(strip(line))
    if lead_char == '('
        vals = parse.(Float64, split(strip(line, ['(', ')']), ','))
        n_count += 1
        return vals[1], vals[2]
    else
        val = parse(Float64, line)
        return val, val
    end
end
include("tjlf_modules.jl")
include("tjlf_geometry.jl")

#
#     parameters: 
#     inputs - InputTJLF struct constructed using the input.TGLF file
#     
#     outputs:
#     ky_spectrum - array of floats that form the ky grid
#     nky - size of the ky_spectrum array
#
#     the input file provides the type of kygrid to create (values 1 to 5)
#     and this function creates it accordingly
# 

function get_ky_spectrum(inputs::InputTJLF{T}, grad_r0::T)::Tuple{Vector{T}, Int} where T<:Real
    
    ### values from input
    units_in = inputs.UNITS
    nky_in = inputs.NKY
    spectrum_type = inputs.KYGRID_MODEL
    rho_e = √(inputs.TAUS[1]*inputs.MASS[1])/ abs(inputs.ZS[1])
    rho_ion = √(inputs.TAUS[2]*inputs.MASS[2])/ abs(inputs.ZS[2])
   
    ### value from modules.f90
    ky_min = 0.05
    ky_max = 0.7
    ky_in = 0.3


    if(units_in == "GYRO")
        ky_factor = 1.0
    else
        ky_factor = grad_r0
    end

    ky0 = ky_min
    ky1 = ky_max
    ### have no idea why this is hard coded like this -DSUN
    nk_zones = 3
    if(nk_zones>=2) ky1 = ky_max/√(inputs.MASS[1]) end


    if(spectrum_type == 0)
        nky = nky_in
        ky_spectrum = Vector{Float64}(undef,nky)
        dky_spectrum = Vector{Float64}(undef,nky)
        ky1 = ky_in
        dky0 = ky1/nky_in
        for i = 1:nky
            ky_spectrum[i] = i*dky0
            dky_spectrum[i] = dky0
        end

    elseif(spectrum_type == 1)   # APS07 spectrum
        nky = 9
        ky_spectrum = Vector{Float64}(undef,nky + nky_in)
        dky_spectrum = Vector{Float64}(undef,nky + nky_in)

        ky_max = 0.9*ky_factor/rho_ion
        dky0 = ky_max/nky
        for i = 1:nky
            ky_spectrum[i] = i*dky0
            dky_spectrum[i] = dky0
        end

        ky0 = ky_max + dky0
        ky1 = 0.4*ky_factor/rho_e

        
        ### will this if ever fail?
        if(nky_in > 0)
            ##### divide by zero?
            dky0 = log(ky1/ky0)/(nky_in-1)
            lnky = log(ky0)
            for i = nky+1:nky+nky_in     
                ky_spectrum[i] = exp(lnky)
                dky_spectrum[i] = ky_spectrum[i]*dky0
                lnky = lnky + dky0
            end
            nky = nky + nky_in
        end

    elseif(spectrum_type == 2)  # IAEA08 spectrum
        nky1 = 8
        nky2 = 7
        nky = nky1 + nky2

        ky_spectrum = Vector{Float64}(undef,nky + nky_in)
        dky_spectrum = Vector{Float64}(undef,nky + nky_in)

        dky0 = 0.05*ky_factor/rho_ion
        for i = 1:nky1
            ky_spectrum[i] = i*dky0
            dky_spectrum[i] = dky0
        end
        dky0 = 0.2/rho_ion
        ky0 = ky_spectrum[nky1]


        for i = nky1+1:nky
            ky_spectrum[i] = ky0 + (i-nky1)*dky0
            dky_spectrum[i] = dky0
        end
        ky0 = ky_spectrum[nky] + dky0
        ky1 = 0.4*ky_factor/rho_e
        
        if(nky_in > 0)
            ###### divide by zero?
            dky0 = log(ky1/ky0)/(nky_in-1)
            lnky = log(ky0)
            for i = nky+1:nky+nky_in     
                ky_spectrum[i] = exp(lnky)
                dky_spectrum[i] = ky_spectrum[i]*dky0
                lnky = lnky + dky0
            end
            nky = nky + nky_in
        end

    elseif(spectrum_type == 3)   # ky_min=ky_in spectrum similar to APS07
        ky_max = 1.0*ky_factor/rho_ion 
        nky1 = Int(floor(ky_max/ky_in)) - 1
        nky2 = 1
        nky = nky1 + nky2
        ky_spectrum = Vector{Float64}(undef,nky + nky_in)
        dky_spectrum = Vector{Float64}(undef,nky + nky_in)

        ky_min = ky_in
        dky0 = ky_min
        ky_spectrum[1] = ky_min
        dky_spectrum[1] = ky_min
        for i = 2:nky1
            ky_spectrum[i] = ky_spectrum[i-1] + dky0
            dky_spectrum[i] = dky0
        end

        if(ky_spectrum[nky1] < ky_max)
            nky2 = 1
            ky_min = ky_spectrum[nky1]
            dky0 = (ky_max-ky_min)/nky2
            ky_spectrum[nky1+1] = ky_spectrum[nky1] + dky0
            dky_spectrum[nky1+1] = dky0
        else
            nky2 = 0
            nky = nky1 + nky2
            ky_max = ky_spectrum[nky1]
            ### truncate the array because it's one less
            ky_spectrum = ky_spectrum[1:nky_in+nky]
        end

        ### gonna add this to make it look the same
        if(nky_in > 0)
            ky0 = ky_max
            ky1 = 0.4*ky_factor/rho_e 
            lnky = log(ky0+dky0)
            for i = nky+1:nky+nky_in
                ky_spectrum[i] = exp(lnky)
                dky_spectrum[i] = ky_spectrum[i]*dky0
                lnky = lnky + dky0
            end
            nky = nky + nky_in
        end

    elseif(spectrum_type == 4)   # APS07 spectrum with ky_min=0.05
        nky1 = 5
        nky2 = 7
        nky = nky1 + nky2

        ky_spectrum = Vector{Float64}(undef,nky + nky_in)
        dky_spectrum = Vector{Float64}(undef,nky + nky_in)

        ky_min = 0.05*ky_factor/rho_ion
        for i = 1:nky1
            ky_spectrum[i] = i*ky_min
            dky_spectrum[i] = ky_min
        end
        ky_min = ky_spectrum[nky1]
        ky_max = 1.0*ky_factor/rho_ion
        
        dky0 = 0.1*ky_factor/rho_ion
        for i = nky1+1:nky
            ky_spectrum[i] = ky_min + (i-4)*dky0
            dky_spectrum[i] = dky0
        end


        if(nky_in > 0)
            ky0 = 1.0*ky_factor/rho_ion
            ky1 = 0.4*ky_factor/rho_e
            dky0 = log(ky1/ky0)/(nky_in-1)
            lnky = log(ky0)
            dky_spectrum[nky]=ky0-ky_spectrum[nky]
            for i = nky+1:nky+nky_in     
                ky_spectrum[i] = exp(lnky)
                dky_spectrum[i] = ky_spectrum[i]*dky0
                lnky = lnky + dky0
            end
            nky = nky + nky_in
        end
    
    elseif(spectrum_type == 5)
        error("apparently no one uses this one, it is also sus -DSUN")
        # nky1 = 10

        # ky1 = ky_in
        # dky0 = ky1/nky_in
        # ntest = Int(floor(0.5/dky0))
        # ####### nky_in can't be too large or else it's bad
        # nky = nky1 + nky_in - ntest
        # if nky<nky1
        #     error("nky_in is too large and IDK what I'm supposed to do,
        #     basically, the ky_spectrum array created will be too small now")
        # end
        # ky_spectrum = zeros(nky)
        # dky_spectrum = zeros(nky)

        # for i = 1:nky1
        #     ky_spectrum[i] = i*0.05
        #     dky_spectrum[i] = 0.05
        # end

        # ky_spectrum[nky1+1] = ntest*dky0
        # #### is this a typo? why are we missing dky_spectrum[nky1+1]
        # dky_spectrum[nky1] = ky_spectrum[nky1+1]-0.5

        # for i=nky1+2:nky
        #     ky_spectrum[i] = ky_spectrum[i-1] + dky0
        #     dky_spectrum[i] = dky0
        # end
    end

    return ky_spectrum, nky

end
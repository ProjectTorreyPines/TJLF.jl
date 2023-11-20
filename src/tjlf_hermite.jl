"""
    function gauss_hermite(inputs::InputTJLF{T}) where T<:Real

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl

outputs:
    OutputHermite(x,wx,h)               - OutputHermite struct

description:
    computes abscisca's and weights for Gauss- Hermite integration adapted from Numerical Recipes in Fortran Pg.147
    reversed order of roots so that x(m) is largest and x(1) is smallest only the positive roots are found and stored
    #--------------------------------------------------------------
    initializes the hermite basis functions and ten point Gauss-Hermite integration with absissas x and weights w.
"""
function gauss_hermite(inputs::InputTJLF{T}) where T<:Real

    eps=3.0E-14
    maxit=100

    # set up the hermite basis x-grid
    nbasis = inputs.NBASIS_MAX
    nx = 2*inputs.NXGRID -1
    h0 = 1.0/π^0.25
    m = Int((nx+1)/2) ### NXGRID

    y = Vector{Float64}(undef, m)
    wy = Vector{Float64}(undef, m)
    z = 0.0
    for i = m:-1:1
        # initial guess z
        if(i==m)
            z = √(2*nx + 1) - 1.85575*(2*nx + 1)^(-0.16667)
        elseif(i==m-1)
            z = z - 1.14*nx^.426/z
        elseif(i==m-2)
            z = 1.86*z - 0.86*y[m]
        elseif(i==m-3)
            z = 1.91*z - 0.91*y[m-1]
        else
            z = 2.0*z - y[i+2]
        end

        pp = 0
        for its = 1:maxit
            p1 = h0
            p2 = 0.0
            for j = 1:nx
                p3 = p2
                p2 = p1
                p1 = z*√(2.0/j)*p2 - √((j-1)/j)*p3
            end
            pp = √(2.0*nx)*p2
            z1 = z
            z = z1 - p1/pp
            if(abs(z-z1)<=eps) break end
        end

        y[i] = z
        wy[i] = 2/pp^2
    end
    wy[1] = wy[1]/2

    x = Vector{Float64}(undef, nx)
    wx = Vector{Float64}(undef, nx)
    for i = 1:m
        x[m+i-1] = y[i]
        x[i] = -y[m+1-i]

        wx[m+i-1] = wy[i]/2
        wx[i] = wy[m+1-i]/2
    end
    x[m] = 0.0
    wx[m] = 2.0*wx[m]

    #--------------------------------------------------------------

    #     The hermite basis functions of the wave function are
    #       H(1,i) = 1.0
    #       H(2,i) = 2.0*x(i)
    #       for j>= 3
    #       H(j,i) = 2.0*x(i)*H(j-1,i)-2.0*(j-2)*H(j-2,i)
    #
    #       H(3,i) = 4.0*x(i)**2 -2.0
    #       H(4,i) = 8.0*x(i)**3 -12.0*x(i)
    #       H(5,i) = 16.0*x(i)**4 -48.0*x(i)**2 + 12.0
    #       H(6,i) = 32.0*x(i)**5 -160.0*x(i)**3 + 120.0*x(i)
    #       H(7,i) = 64.0*x(i)**6 -480.0*x(i)**4 +720.0*x(i)**2 -120.0
    #       H(8,i) = 128.0*x(i)**7 -1344.0*x(i)**5 +3360.0*x(i)**3 -1680.0*x(i)
    #       H(9,i) = 256.0*x(i)**8 -3584.0*x(i)**6 +13440.0*x(i)**4 -13440.0*x(i)**2 +1680.0
    #       H(10,i)= 512.0*x(i)**9 -9216.0*x(i)**7 +48384.0*x(i)**5 -80640.0*x(i)**3 +30240.0*x(i)
    #       H(11,i)= 1024.0*x(i)**10 -23040.0*x(i)**8 +161280.0*x(i)**6 -403200.0*x(i)**4 +302400.0*x(i)**2 -30240.0
    #     The normalized hermite polynomial are
    #       h(j+1,i) = H(j,i)/SQRT(sum(w*H(j)**2))
    #     These have been pre-evaluated for speed.
    #     They are orthonormal with error of 1.0E-8.

    h0 = √(2)/π^0.25
    h = Matrix{Float64}(undef, nbasis, nx)
    h[1,:] .= h0
    h[2,:] .= x.*(√(2)*h0)

    if nbasis>2
        for i = 1:nx
            for j = 3:nbasis
                h[j,i] = x[i]* √(2.0/(j-1)) * h[j-1,i] - √((j-2)/(j-1)) * h[j-2,i]
            end
        end
    end

    return OutputHermite(x,wx,h,size(inputs.KY_SPECTRUM,1))
end
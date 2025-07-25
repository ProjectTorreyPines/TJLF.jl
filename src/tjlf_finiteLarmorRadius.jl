const nf = 40
const na = 13
# const h = Vector{Float64}(undef, na) #### CAUSES ERROR IN MULTITHREADING

function FLR_hs(ft::T, b::T, Amat::Matrix{T}; b2=b^2, b25=b^2.5)::T where T<:Real
#******************************************************************
#     Approxmation to the integral of (J0^2)*Fmaxwellian over a
#     wedge of velocity space -ft < (v_par/v) < ft.
#     b = (k_per*√(T/m)/(eB/mc))^2
#     The approximation has the form
#     Hn = ft*(h0+a1*h1+a2*h2+a3*h3+a4*h4+a5*h5+a6*h6
#             +a7*h7+a8*h8+a9*h9+a10*h10+a11*h11++a12*h12+a13*h13)
#     where a1 - a13 are fit coefficients which are
#     functions of ft and
#     xa1 = (0.25)^4
#     xa2 = (0.5)^4
#     xa3 = (0.75)^4
#     xa4 = (1.0)^4
#     xa5 = 0.25*(1.5)^5
#     xa6 = 0.25*(2.0)^5
#     xa7 = 0.25*(2.5)^5
#     xa8 = 0.25*(3.0)^5
#     xa9 = 0.25*(4.0)^5
#     xa10 = 0.25*(6.0)^5
#     xa11 = 0.25*(9.0)^5
#     xa12 = 0.25*(15.0)^5
#     xa13 = 0.25*(24.0)^5
#     h0 = 1/(1+b)
#     h1 = b/(xa1+b^2)
#     h2 = b/(xa2+b^2)
#     h3 = b/(xa3+b^2)
#     h4 = b/(xa4+b^2)
#     h5 = b^2/(xa5+b^2.5)
#     h6 = b^2/(xa6+b^2.5)
#     h7 = b^2/(xa7+b^2.5)
#     h8 = b^2/(xa8+b^2.5)
#     h9 = b^2/(xa9+b^2.5)
#     h10 = b^2/(xa10+b^2.5)
#     h11 = b^2/(xa11+b^2.5)
#     h12 = b^2/(xa12+b^2.5)
#     h13 = b^2/(xa13+b^2.5)
#******************************************************************

    y4 = FLR_constants.y4
    y5 = FLR_constants.y5
    g = FLR_constants.g

    h14 = i -> b/(y4[i] + b2)
    h5na = i -> b2/(0.25*y5[i] + b25)

    # transform to gt grid
    gt = √(1-ft)
    if(gt>g[nf]) gt = g[nf] end
    # find grid position
    if(gt<=g[2])
        i = 1
    else
        i = round(Int, gt/0.025, RoundDown)
    end
    j = i+1

    # interpolate the coefficients
    dg = (gt - g[i])/(g[j]-g[i])
    hs = 0.
    for k = 1:na
        ca = Amat[i,k] + (Amat[j,k] - Amat[i,k])*dg
        # sum up the terms
        hs = hs + ca*(k <= 4 ? h14(k) : h5na(k))
    end

    return hs
end

#------------------------------------------------------------------
function FLR_Hn(ft::T, b::T; b2=b^2, b25=b^2.5)::T where T<:Real
    hs = 1.0/(1.0+b) + FLR_hs(ft, b, FLR_constants.a_Hn)
    return ft*hs
end

#------------------------------------------------------------------
function FLR_dHp1(ft::T, b::T; b2=b^2, b25=b^2.5)::T where T<:Real
    return (ft^3)*FLR_hs(ft, b, FLR_constants.a_dHp1)
end

#------------------------------------------------------------------
function FLR_dHp3(ft::T, b::T; b2=b^2, b25=b^2.5)::T where T<:Real
   return ft*FLR_hs(ft, b, FLR_constants.a_dHp3)
end

#------------------------------------------------------------------
function FLR_dHr11(ft::T, b::T; b2=b^2, b25=b^2.5)::T where T<:Real
    return (3.0*ft^5)*FLR_hs(ft, b, FLR_constants.a_dHr11)
end

#------------------------------------------------------------------
function FLR_dHr13(ft::T, b::T; b2=b^2, b25=b^2.5)::T where T<:Real
    return (5.0/3.0)*(ft^3)*FLR_hs(ft, b, FLR_constants.a_dHr13)
end

#------------------------------------------------------------------
function FLR_dHr33(ft::T, b::T; b2=b^2, b25=b^2.5)::T where T<:Real
    return (5.0/3.0)*ft*FLR_hs(ft, b, FLR_constants.a_dHr33)
end

#------------------------------------------------------------------
function FLR_dHw113(ft::T, b::T; b2=b^2, b25=b^2.5)::T where T<:Real
    return 7.0*ft^5*FLR_hs(ft, b, FLR_constants.a_dHw113)
end

#------------------------------------------------------------------
function FLR_dHw133(ft::T, b::T; b2=b^2, b25=b^2.5)::T where T<:Real
    return (35.0/9.0)*ft^3*FLR_hs(ft, b, FLR_constants.a_dHw133)
end

#------------------------------------------------------------------
function FLR_dHw333(ft::T, b::T; b2=b^2, b25=b^2.5)::T where T<:Real
      return (35.0/9.0)*ft*FLR_hs(ft, b, FLR_constants.a_dHw333)
end



###### one dot . is within the module, two dots .. is for outside the moduel



module temp
include("test2.jl")
using .test2
println(nb)
global test2.nb += 3
println(nb)
end

A = [1,2,3]
println(A[3])



module ModuleA
    export var

    const global var = "Hello from ModuleA!"
end

module ModuleB
    # import ..ModuleA: var  # Import the variable var from ModuleA using a double dot
    using ..ModuleA
    println(ModuleA.var)
    global ModuleA.var = "Hello from ModuleB!" # gives a warning
end

function test()
    println(ModuleA.var) # fine
    # println(var) # errror
end

test()

function scope()
    # a = [1,2,3,4]
    # b = [2,3,4,5]
    # c = zeros(4)

    # c .= ifelse.(b.>a, b.+max.(b,a.*3), 0)
    # println(c)

    if false
        a = "asdf"
    else
        a = 2
    end
    print(a)

end

scope()

matrixZ = Vector{ComplexF64}(undef, 120*120)
matrixV = Vector{ComplexF64}(undef, 120)

zmatDirectory = "../../zmat"
zmat = readlines(zmatDirectory)
vDirectory = "../../v"
v = readlines(vDirectory)
index = 1
for i in eachindex(zmat)
    if contains(zmat[i], "=") continue end
    zLine = split(strip(zmat[i], ['(', ')',' ']),",")
    zLineR = parse(Float64, zLine[1])
    zLineI = parse(Float64, zLine[2])
    matrixZ[index] = zLineR + im*zLineI
    index+=1
end
matrixZ = Matrix(reshape(matrixZ,(120,120))')
index = 1
for i in eachindex(v)
    if contains(v[i], "=") continue end
    vLine = split(strip(v[i], ['(', ')',' ']),",")
    vLineR = parse(Float64, vLine[1])
    vLineI = parse(Float64, vLine[2])
    matrixV[index] = vLineR + im*vLineI
    index+=1
end
v = zeros(ComplexF64, 120)
v .= 1.0e-13
gesv!(matrixZ,v)
v = matrixZ \ v
norm(v)^2
cond(matrixZ)
matrixZ * v
matrixZ * matrixV




zmatDirectory = "../../zmat"
zmat = readlines(zmatDirectory)
zmatDirectory2 = "../../zmat2"
zmat2 = readlines(zmatDirectory2)
for i in eachindex(zmat)
    if contains(zmat[i], "=") continue end
    zLine = split(strip(zmat[i], ['(', ')',' ']),",")
    zLineR = parse(Float64, zLine[1])
    zLineI = parse(Float64, zLine[2])

    try
        zLine2 = parse(ComplexF64, strip(zmat2[i], [' ','[',']']))
        zLineR2 = real(zLine2)
        zLineI2 = imag(zLine2)
        if !isapprox(zLineR,zLineR2,atol=10^-12)
            println(i)
            println(zLineR)
            println(zLineR2)
        end
        if !isapprox(zLineI,zLineI2,atol=10^-12)
            println(i)
            println(zLineI)
            println(zLineI2)
        end
    catch
        println(i)
        error("FAIL")
    end


end
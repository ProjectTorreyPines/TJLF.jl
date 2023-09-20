function a()


asdf=0
println(asdf)
println(2π
+3)
if(true) println("one line") end
end

module test2
export nb
nb = 1

end

arr=[1,2,3,4,5,6]
b=[1,2,3,4,5,6]
# c=zero(6)
arr[1]+=1
c =arr.*b
c[1] = 4
# println(c)
# println(arr)
# println(.√(b.^2 .+ b))

d = zeros(6)
d .= arr
d[1]+=1
# println(arr)

e = 2 ./b
# println(e)

f = [1,2,3,4]
# f .= ifelse.(f.<=2,f, f.-1)
println(ifelse.(f.<=2,f, f.-1))

a()













libpath = "./src/simplemodule.so"

h = Int32[3]
ccall((:main, "./src/simplemodule.so"), Cvoid, (Ptr{Int32},), h)

# ccall((:__simplemodule_MOD_foo, lib), Int32, (Ptr{Int32},), h)
# ccall((:__tglf, "./src/Fortran/tglf.so"), Int32, (Ptr{Int32},), h)






# h = "../test_SAT3/"
# ccall((:__tglf_MOD_foo, "./src/Fortran/tglf.so"), Cvoid, (Ptr{Cchar},), h)



ccall((:main, "./src/Fortran/tglf.so"), Cvoid, () ,)
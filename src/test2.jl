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
println(c)
println(arr)
println(.√(b.^2 .+ b))

d = zeros(6)
d .= arr
d[1]+=1
println(arr)

a()
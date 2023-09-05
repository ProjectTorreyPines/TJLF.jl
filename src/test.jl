

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
#=========================================================================#
# The following portion of the main.jl is mean for translating driver.f90 #
#=========================================================================#
using MPI
#function driver()
    comm = MPI.COMM_WORLD
    # There are some fundamental parts of Distributed that I could use. I could also try to use MPI.jl?
    MPI.Init()
    MPI.Comm_rank(comm)
    MPI.Comm_size(comm) # I have only one thing I can access here, so this isn't shocking.

    println("Hello world, I am $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm)) \n")
#end
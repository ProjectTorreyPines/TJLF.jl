
A = rand(10, 10) * 10

using LinearAlgebra

A = LinearAlgebra.Symmetric(A)

for i in 1:10

    # using CPU
    #############################################################################

    eig = @time eigen(A)
    lambdas = eig.values
    evecs = eig.vectors

    vec1 = evecs[:, 1]
    vec2 = evecs[:, 2]

    lhs1 = A * vec1
    rhs1 = lambdas[1] * vec1

    # println("Eigenvalues: $lambdas")
    # println("Eigenvectors: $evecs")

    # using GPU
    #############################################################################

    using CUDA

    A_gpu = CuArray(A)
    eig_gpu = @time CUSOLVER.syevd!('V', 'U', A_gpu)
    lambdas_gpu = eig_gpu[1]
    evecs_gpu = eig_gpu[2]

    vec1_gpu = evecs_gpu[:, 1]
    vec2_gpu = evecs_gpu[:, 2]

    lhs1_gpu = A_gpu * vec1_gpu
    rhs1_gpu = lambdas_gpu[1:1] .* vec1_gpu

    # println("GPU Eigenvalues: $lambdas_gpu")
    # println("GPU Eigenvectors: $evecs_gpu")

end

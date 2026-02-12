
A = rand(1000, 1000) * 10

using LinearAlgebra
@time import CUDA

# Preallocate arrays to reuse (avoid allocations in loop)
A_cpu = copy(A)  # Copy for CPU (syevd! is destructive)
A_gpu = CUDA.CuArray(A)  # Move to GPU once

# Preallocate output arrays
lambdas = Vector{Float64}(undef, size(A, 1))
evecs = Matrix{Float64}(undef, size(A))
lambdas_gpu = CUDA.CuArray{Float64}(undef, size(A, 1))
evecs_gpu = CUDA.CuArray{Float64}(undef, size(A))

for i in 1:10

    # using CPU
    #############################################################################
    
    copyto!(A_cpu, A)  # Reset A_cpu (syevd! modifies it)
    eig = @time LAPACK.syevd!('V','U',A_cpu)

    # Use views to avoid allocations
    # vec1 = @view eig[2][:, 1]
    # vec2 = @view eig[2][:, 2]

    # lhs1 = A * vec1
    # rhs1 = eig[1][1] * vec1

    # println("Eigenvalues: $lambdas")
    # println("Eigenvectors: $evecs")

    # using GPU
    #############################################################################

    copyto!(A_gpu, A)  # Reset A_gpu (syevd! modifies it)
    
    # Synchronize before timing for accurate measurement
    CUDA.synchronize()
    t_start = time()
    
    eig_gpu = @time CUDA.CUSOLVER.syevd!('V', 'U', A_gpu)
    
    # Synchronize after to ensure operation completes
    CUDA.synchronize()
    t_end = time()
    
    # println("  $(t_end - t_start) seconds (GPU synchronized)")
    
    # Use views to avoid allocations
    # vec1_gpu = @view eig_gpu[2][:, 1]
    # vec2_gpu = @view eig_gpu[2][:, 2]

    # lhs1_gpu = A_gpu * vec1_gpu
    # rhs1_gpu = eig_gpu[1][1:1] .* vec1_gpu

    # println("GPU Eigenvalues: $lambdas_gpu")
    # println("GPU Eigenvectors: $evecs_gpu")

    # add @view to decrease number of allocations

end

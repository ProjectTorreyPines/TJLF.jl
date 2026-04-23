__precompile__(false)
module TJLFCUDAExt

using TJLF
using CUDA

TJLF._cuda_functional() = CUDA.functional()
TJLF._cuda_device_count() = length(CUDA.devices())

function TJLF._gpu_solve!(A::Matrix{ComplexF64}, B::Matrix{ComplexF64})
    A_gpu = CUDA.CuArray(A)
    B_gpu = CUDA.CuArray(B)
    ipiv = CUDA.CuArray{Int32}(undef, size(B_gpu, 1))
    B_factored, ipiv, _ = CUDA.CUSOLVER.getrf!(B_gpu, ipiv)
    A_gpu = CUDA.CUSOLVER.getrs!('N', B_factored, ipiv, A_gpu)
    return Array(A_gpu)
end

end

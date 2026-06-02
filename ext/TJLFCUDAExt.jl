module TJLFCUDAExt

using CUDA
using TJLF

# Per-device dispatch locks. cusolverDnCreate races when many threads call it
# simultaneously on first use, so Julia-level dispatch is serialized — but only
# *per device*. Using one lock per GPU (rather than a single global lock) lets
# solves on different GPUs proceed concurrently, which is what enables
# multi-GPU-per-radius scaling: e.g. with --gpus-per-task=2 the round-robin below
# splits a radius's ~1024 eigensolves across both A100s, ~halving per-radius time.
# Within a single device, dispatch stays serialized (the GPU is internally
# parallel, so this has negligible cost). With one visible GPU this reduces to a
# single lock — identical to the previous behavior.
const _LOCKS_GUARD = ReentrantLock()
const _GPU_DEVICE_LOCKS = Dict{Int,ReentrantLock}()
function _device_lock(id::Integer)
    lock(_LOCKS_GUARD) do
        get!(() -> ReentrantLock(), _GPU_DEVICE_LOCKS, Int(id))
    end
end

# Round-robin counter + cached device list to spread successive solves across all
# visible GPUs. Lazily initialised (double-checked) so it works whether or not
# CUDA is functional at load time.
const _GPU_RR = Threads.Atomic{Int}(0)
const _DEVICES = Ref{Vector{CuDevice}}(CuDevice[])
function _devices()
    devs = _DEVICES[]
    isempty(devs) || return devs
    lock(_LOCKS_GUARD) do
        devs = _DEVICES[]
        isempty(devs) || return devs
        devs = collect(CUDA.devices())
        _DEVICES[] = devs
        return devs
    end
end

function __init__()
    TJLF._CUDA_FUNCTIONAL[]   = () -> CUDA.functional()
    TJLF._CUDA_DEVICE_COUNT[] = () -> length(CUDA.devices())
    # Returns Vector{ComplexF64} of eigenvalues of B⁻¹A, computed entirely on GPU
    # (CUSOLVER getrf/getrs for the solve, then Xgeev for the diagonalisation).
    # B is only needed for the LU factorisation; its overwritten form is unused.
    TJLF._CUDA_SOLVE[] = (A::Matrix{ComplexF64}, B::Matrix{ComplexF64}) -> begin
        devs = _devices()
        n = length(devs)
        # round-robin device selection across all visible GPUs
        dev = devs[mod1(Threads.atomic_add!(_GPU_RR, 1) + 1, n)]
        CUDA.device!(dev)
        lock(_device_lock(CUDA.deviceid(dev))) do
            A_gpu = CUDA.CuArray(A)
            B_gpu = CUDA.CuArray(B)
            ipiv  = CUDA.CuArray{Int32}(undef, size(B_gpu, 1))
            B_factored, ipiv, _ = CUDA.CUSOLVER.getrf!(B_gpu, ipiv)
            A_gpu = CUDA.CUSOLVER.getrs!('N', B_factored, ipiv, A_gpu)
            W, _, _ = CUDA.CUSOLVER.Xgeev!('N', 'N', A_gpu)
            return Array(W)
        end
    end
end

end

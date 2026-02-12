
using JLD2

# CPU = load("run_tjlf_CPU.jld2")
# GPU = load("run_tjlf_GPU.jld2")
# CPU_GPU_node = load("run_tjlf_CPU(GPU_node).jld2")

# Cavg = CPU["avg"]
# Gavg = GPU["avg"]
# CGavg = CPU_GPU_node["avg"]

# Ctimes = CPU["times"]
# Gtimes = GPU["times"]
# CGtimes = CPU_GPU_node["times"]

# Cwarmup = CPU["time_warmup"]
# Gwarmup = GPU["time_warmup"]
# CGwarmup = CPU_GPU_node["time_warmup"]

using Plots

G_avgs = Array{Float64}(undef, 15)
CG_avgs = Array{Float64}(undef, 15)

nbasis = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 28, 30]

for n in nbasis
    # Cdata = load("outputs/CPU/run_tjlf_CPU_$n.jld2")
    Gdata = load("outputs/GPU/run_tjlf_GPU_$n.jld2")
    CGdata = load("outputs/CPU(GPU_node)/run_tjlf_CPU(GPU_node)_$n.jld2")

    # Cavg = Cdata["avg"]
    Gavg = Gdata["avg"]
    CGavg = CGdata["avg"]

    # Cwarmup = Cdata["time_warmup"]
    # Gwarmup = Gdata["time_warmup"]
    # CGwarmup = CGdata["time_warmup"]

    percent_imrpovement = 100 - (Gavg / CGavg * 100)
    G_avgs[findfirst(isequal(n), nbasis)] = Gavg
    CG_avgs[findfirst(isequal(n), nbasis)] = CGavg

    println("NBASIS = $n: GPU avg: $Gavg s, CPU (GPU node) avg: $CGavg s, percent improvement: $percent_imrpovement %")
    # println("NBASIS = $n: GPU warmup: $Gwarmup s, CPU (GPU node) warmup: $CGwarmup s")
end

plot(nbasis, G_avgs, label="GPU")
plot!(nbasis, CG_avgs, label="CPU (GPU node)")

savefig("outputs/GPU_CGPU_vs_NBASIS.png")

# println("avgs: CPU: $Cavg s, GPU: $Gavg s, CPU (GPU node): $CGavg s")
# println("warmups: CPU: $Cwarmup s, GPU: $Gwarmup s, CPU (GPU node): $CGwarmup s")

include("../src/MCWF.jl")
include("../src/model.jl")

# パラメータの設定
params = Parameters(Δ=0.0, Γ=1/6, Ω=1.0)

H=makeHamiltonian(params)
L=makeLindbladian(params)

ψ0=zeros(ComplexF64, 2)
ψ0[2]=1.0

dt=0.01
tmax=20

ts,ψs1 = trajectory_approach.time_evol_trajectory(ψ0, dt, tmax, H, [L])
ts,ψs2 = trajectory_approach.time_evol_trajectory(ψ0, dt, tmax, H, [L])

using Plots, LaTeXStrings

fig=plot(ts, [abs(ψs1[i][1])^2 for i ∈ 1:length(ts)], xlabel=L"\Omega t", ylabel=L"P_{\mathrm{e}}", label="MCWF (sample 1)", ylims=(0,1), lw=2)
plot!(ts, [abs(ψs2[i][1])^2 for i ∈ 1:length(ts)], label="MCWF (sample 2)", ylims=(0,1), lw=2)
hline!(fig, [steadyPe(params)], label="", ls=:dash)
savefig(fig, "fig2_MCWF_samples.pdf")
savefig(fig, "fig2_MCWF_samples.png")

# 100サンプルの実行と平均の計算
n_samples = 100
tmax=40
ts=0:dt:tmax
all_Pe = zeros(Float64, length(ts), n_samples)

for i in 1:n_samples
    _, ψs = trajectory_approach.time_evol_trajectory(ψ0, dt, tmax, H, [L])
    all_Pe[:,i] = [abs(ψs[j][1])^2 for j in 1:length(ts)]
end

using Statistics
# 平均値の計算
avg_Pe = [mean(all_Pe[i,:]) for i in 1:length(ts)]
# 標準偏差の計算
std_Pe = [std(all_Pe[i,:])/sqrt(n_samples) for i in 1:length(ts)]

# 平均値のプロット
fig2=plot(ts, avg_Pe, ribbon=std_Pe, fillalpha=0.5, label="MCWF (average of $n_samples samples)", lw=2)
hline!(fig2, [steadyPe(params)], label="", ls=:dash)
savefig(fig2, "fig3_MCWF_average.pdf")
savefig(fig2, "fig3_MCWF_average.png")



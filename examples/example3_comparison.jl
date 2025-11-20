include("../src/ED.jl")
include("../src/MCWF.jl")
include("../src/model.jl")

using Main.sparse_Liouville_space, Statistics, Plots, LaTeXStrings

# パラメータの設定
params = Parameters(Δ=0.0, Γ=1/6, Ω=1.0)

H = makeHamiltonian(params)
L = makeLindbladian(params)

# 初期状態の設定
ρ0 = zeros(ComplexF64, 2, 2)
ρ0[2,2] = 1.0
ψ0 = zeros(ComplexF64, 2)
ψ0[2] = 1.0

# 時間発展のパラメータ
dt = 0.01
tmax = 40
ts = 0:dt:tmax

# EDによる時間発展
ρs = sparse_Liouville_space.time_evol_ED(H, [1.0], [L], ρ0, ts)
Pe_ED = [real(ρ[1,1]) for ρ in ρs]

# MCWFによる時間発展（1000サンプル）
n_samples = 1000
all_Pe = zeros(Float64, length(ts), n_samples)

for i in 1:n_samples
    _, ψs = trajectory_approach.time_evol_trajectory(ψ0, dt, tmax, H, [L])
    all_Pe[:,i] = [abs(ψs[j][1])^2 for j in 1:length(ts)]
end

# MCWFの平均値と標準誤差の計算
avg_Pe = [mean(all_Pe[i,:]) for i in 1:length(ts)]
std_Pe = [std(all_Pe[i,:])/sqrt(n_samples) for i in 1:length(ts)]

# プロット
fig = plot(ts, avg_Pe, ribbon=std_Pe, fillalpha=0.5, lw=2, label="MCWF ($n_samples samples)")
plot!(ts, Pe_ED, label="ED", lw=2, ls=:dash, xlabel=L"\Omega t", ylabel=L"P_{\mathrm{e}}")
hline!([steadyPe(params)], label="steady state", ls=:dash)
savefig(fig, "../assets/fig4_ED_vs_MCWF.pdf")
savefig(fig, "../assets/fig4_ED_vs_MCWF.png")

include("ED.jl")
include("model.jl")

using Main.sparse_Liouville_space

# パラメータの設定
params = Parameters(Δ=0.0, Γ=1/6, Ω=1.0)

# 初期状態の設定（基底状態）
ρ0 = sparse([0.0 0.0; 0.0 1.0])

# 時間の設定
ts = range(0, 40, length=1000)

# ハミルトニアンとリンドブラディアンの生成
H = makeHamiltonian(params)
L = makeLindbladian(params)

# 時間発展の計算
ρs = sparse_Liouville_space.time_evol_ED(H, [1.], [L], ρ0, ts)

using Plots, LaTeXStrings

fig=plot(ts, real.([ρs[i][1,1] for i ∈ 1:length(ts)]), xlabel=L"\Omega t", ylabel=L"\rho_{11}", label="ED", ylims=(0,1), lw=2)
hline!(fig, [steadyPe(params)], label="", ls=:dash)
savefig(fig, "fig1_ED.pdf")
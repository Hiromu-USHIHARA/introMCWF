この記事では，[前回](https://zenn.dev/hiromu_ushihara/articles/db67c435b9b95b)に続いて，開放量子系の数値計算手法について扱います．
特に今回はMonte Carlo Wave Function (MCWF)法の導入・実装と厳密対角化（ED）との比較を行います．

> このチュートリアルは[Zennでも公開](https://zenn.dev/hiromu_ushihara/articles/c6ef07f16666ee)していますので，数式がうまく表示されない場合はそちらをご覧ください．



## 扱う問題（再確認）

このチュートリアルでは，Gorini--Kosakkowski--Sudarshan--Lindblad (GKSL)方程式[1--3]

$$
\frac{d\rho}{dt} = -i \left[ H, \rho(t) \right]_{-} + \sum_i \left( L_i \rho L_i^\dagger - \frac{1}{2} \left[ L_i^\dagger L_i, \rho \right]_{+} \right).
$$

で記述される開放量子系を扱います．
特に具体例として，以下で定義される電磁場と結合した二準位系を扱います[5]:

$$H=-\dfrac{\Omega}{2}\sigma_{\mathrm x}-\Delta\sigma_+\sigma_-,$$
$$L=\sqrt\Gamma\sigma_-.$$

## 厳密対角化の問題点

前回の記事で扱った厳密対角化（ED）はストラテジーが明瞭で，かつ数値対角化の精度の範囲で厳密な数値計算手法でした．
しかし，数値計算のコストという観点では問題があります．

この問題は以下のようにして理解出来ます．
具体例として，$N$個のスピン$1/2$からなる系を考えます．
$1$サイトあたりのHilbert空間の次元は$2$で，$N$サイト全体で$2^N$になります．
密度行列のサイズは$2^{N}\times2^N$です．
対角化の対象となる超演算子は演算子を演算子にマップするものでしたから，そのサイズは$2^{2N}\times 2^{2N}$となります．
複素数の要素$1$つあたり$8\times2$バイト必要として，例えば$N=10$とすると，

$$8\times2\times2^{20}\times2^{20}=16\times(2^{10})^4\gtrsim10^{13},$$

すなわち行列を保持するために$10\,\mathrm{TB}$のメモリが必要になります．

もちろん，実際には疎行列を活用するなどの効率化は可能ですが，システムサイズについて指数的にコストが増大するという問題は無視できません．
このことから，**EDでは大規模な系の計算は難しい**ということが分かります．

## モンテカルロ波動関数法（MCWF）

そこで有効となる方法（の一つ）がMCWF [4,5]です^[この記事では一次精度のMCWFのみを扱います．高次精度の方法について知りたい方は文献[5]をご参照ください．]．

### MCWFの手順

まず，GKSL方程式を次のように書き換えます:

$$
\frac{d\rho}{dt} = -i \left(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^\dag\right) + \sum_i L_i \rho L_i^\dagger.
$$

ここで，nonhermitian有効Hamiltonian $H_{\mathrm{eff}}$を導入しました:

$$H_{\mathrm{eff}}:=H-\dfrac{1}{2}\sum_iL_i^\dag L_i.$$

また，初期状態$\left|\phi(t=0)\right>$を純粋状態として用意します.

すると，$\left|\phi(t)\right>$から$\left|\phi(t+\delta t)\right>$への時間発展は次の手順で決定されます．

1. 有効Hamiltonian $H_{\mathrm{eff}}$による時間発展を計算します:
    $$\left|\phi^{(1)}(t+\delta t)\right>=\left(1-iH_{\mathrm{eff}}\right)\left|\phi(t)\right>$$
   > この状態のノルムは
   > $$\left<\phi^{(1)}(t+\delta t)\middle|\phi^{(1)}(t+\delta t)\right>=1-\delta p+O\left(\delta t^2\right)=1-\sum_i\delta p_i+O\left(\delta t^2\right)=1-\delta t\sum_i\left\|L_i\left|\phi(t)\right>\right\|^2+O\left(\delta t^2\right)$$
   > となります．
2. 乱数を生成して，以下の確率的な処理を行います:
   - 確率$1-\delta p$で
     $$\left|\phi(t+\delta t)\right>=\dfrac{\left|\phi^{(1)}(t+\delta t)\right>}{\sqrt{1-\delta p}}$$
   - 確率$\delta p_i/$で
     $$\left|\phi(t+\delta t)\right>=\dfrac{L_i\left|\phi(t)\right>}{\sqrt{\delta p_i/\delta t}}$$
    > 上のケースをnonhermitian time evolution，下のケースをjumpなどと呼ぶことがあります．

これによって規格化された状態$\left|\phi(t+\delta t)\right>$が得られるので，小さな時間刻みに対してこの手順を繰り返すことで時間発展の計算が出来ます．

### MCWFの統計性



### Juliaでの実装

`src/MCWF.jl`

```julia
module trajectory_approach
    using LinearAlgebra, SparseArrays, Random
    # %%
    function time_evol_trajectory(ψ0, dt, tmax, H, Ls)
        Heff = H - 0.5*im*sum(L'*L for L in Ls)
        # es, vs = eigen(Matrix(Heff))
        ts = 0:dt:tmax
        ψs = [Vector{ComplexF64}(ψ0)]
        for i ∈ 1:(length(ts)-1)
            ψ = ψs[i]
            # Vinvψ =  vs\ψ
            # ψ = sum(exp(-im*dt*es[i])*Vinvψ[i]*vs[:,i] for i ∈ 1:size(vs,2))
            ψ = (I(size(H,1)) - im*dt*Heff)*ψ
            δps = [dt*real(ψs[i]'*Ls[j]'*Ls[j]*ψs[i]) for j ∈ 1:length(Ls)]
            δp = sum(δps)
            r1 = rand()
            if r1 > δp
                # ψ = ψ/√(1-δp)
                ψ = ψ/√(ψ'*ψ)
            else
                r2 = rand()
                j = findfirst(cumsum(δps/δp) .> r2)
                ψ = Ls[j]*ψs[i]
                ψ = ψ/√(ψ'*ψ)
            end
            append!(ψs, [ψ])
        end
        return ts, ψs
    end  
end
```

`examples/example2_MCWF.jl`

```julia
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
```
![](https://storage.googleapis.com/zenn-user-upload/f2a2539f23da-20250529.png)


```
# examples/example2_MCWF.jl
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
```

![](https://storage.googleapis.com/zenn-user-upload/36dd0035db02-20250529.png)

## EDとMCWFの比較

`examples/example3_comparison.jl`

```julia
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

# MCWFによる時間発展（100サンプル）
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
```

![](https://storage.googleapis.com/zenn-user-upload/45b058005a2c-20250529.png)



## 参考文献
1. G. Lindblad, "On the generators of quantum dynamical semigroups", Commun. Math. Phys. 48, 119-130 (1976).
2. V. Gorini, A. Kossakowski, and E. C. G. Sudarshan, "Completely positive dynamical semigroups of N-level systems", J. Math. Phys. 17, 821 (1976).
3. H.-P. Breuer, and F. Petruccione, "The Theory of Open Quantum Systems" (Oxford University Press, 2002).
4. J. Dailbard, Y. Castin, and K. Mølmer, "Wave-function approach to dissipative processes in quantum optics", Phys. Rev. Lett. 68, 580 (1992).
5. A. J. Daley, "Quantum trajectories and open many-body quantum systems", Adv. in Phys. 63, 77-149 (2014).
6. J. A. Gyamfi, "Fundamentals of quantum mechanics in Liouville space", Eur. J. Phys. 41, 063002 (2020).

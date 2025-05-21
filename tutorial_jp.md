# 開放量子系の数値計算入門

このチュートリアルでは，開放量子系の数値シミュレーション手法である厳密対角化（exact diagonalization; ED）とモンテカルロ波動関数法（Monte Carlo wave function method; MCWF）を導入・実装します．

## 理論的背景

このチュートリアルで扱うのは開放量子系と呼ばれるタイプの物理系です．開放量子系は着目する系が外部の系に接触しているものを指します．
例えば，熱浴と接触している物質や電磁場と相互作用する原子は開放量子系の例です．

開放量子系の記述として広く用いられているものに，Gorini--Kossakowski--Sudarshan--Lindblad方程式（GKSL方程式）があります^[適当な仮定のもとで開放量子系のダイナミクスがGKSL方程式に従うことが示されます．]:

$$
\frac{d\rho}{dt} = -i\left[H, \rho\right]_- + \sum_i \left( L_i \rho L_i^\dagger - \frac{1}{2}\left[L_i^\dagger L_i, \rho\right]_+ \right)
$$



> 注: このチュートリアルはGPTを使用して作成した下書きです。

## はじめに

このチュートリアルでは、量子系の時間発展をシミュレーションする2つの方法、Exact Diagonalization（ED）法とMonte Carlo Wave Function（MCWF）法を実装します。特に、散逸のある量子系（open quantum systems）に焦点を当て、両手法の特徴と実装方法を解説します。

## 理論的背景

### 散逸のある量子系

散逸のある量子系では、系と環境の相互作用により、量子状態の時間発展は以下のリンドブラッド方程式で記述されます：

\[
\frac{d\rho}{dt} = -i[H, \rho] + \sum_i \left( L_i \rho L_i^\dagger - \frac{1}{2}\{L_i^\dagger L_i, \rho\} \right)
\]

ここで、\(H\)はハミルトニアン、\(L_i\)はリンドブラッド演算子、\(\rho\)は密度行列です。

### ED法（Exact Diagonalization）

ED法は、リンドブラッド方程式を直接数値的に解く方法です。密度行列をベクトル化し、スーパーオペレータを対角化することで時間発展を計算します。

### MCWF法（Monte Carlo Wave Function）

MCWF法は、量子軌道（quantum trajectory）のアンサンブル平均として密度行列を計算する確率的手法です。各軌道は以下のシュレーディンガー型方程式に従います：

\[
\frac{d|\psi\rangle}{dt} = -iH_{\text{eff}}|\psi\rangle
\]

ここで、\(H_{\text{eff}} = H - \frac{i}{2}\sum_i L_i^\dagger L_i\)は有効ハミルトニアンです。

## 実装

### 1. モデルの定義（model.jl）

まず、2準位系のモデルを定義します：

```julia
struct Parameters
    Δ::Float64  # デチューニング
    Γ::Float64  # 減衰率
    Ω::Float64  # ラビ周波数
end

# ハミルトニアンの生成
function makeHamiltonian(params::Parameters)
    H = zeros(ComplexF64, 2, 2)
    H[1,2] = params.Ω/2
    H[2,1] = params.Ω/2
    H[2,2] = params.Δ
    return H
end

# リンドブラッド演算子の生成
function makeLindbladian(params::Parameters)
    L = zeros(ComplexF64, 2, 2)
    L[1,2] = sqrt(params.Γ)
    return L
end
```

### 2. ED法の実装（ED.jl）

ED法では、密度行列をベクトル化し、スーパーオペレータを構築します：

```julia
function time_evol_ED(H, γs, Ls, ρ0, ts)
    # スーパーオペレータの構築
    L = build_superoperator(H, γs, Ls)
    
    # 初期状態のベクトル化
    ρ0_vec = vec(ρ0)
    
    # 時間発展の計算
    ρs = []
    for t in ts
        ρ = exp(L*t) * ρ0_vec
        push!(ρs, reshape(ρ, size(ρ0)))
    end
    return ρs
end
```

### 3. MCWF法の実装（MCWF.jl）

MCWF法では、量子軌道の時間発展を計算します：

```julia
function time_evol_trajectory(ψ0, dt, tmax, H, Ls)
    # 有効ハミルトニアンの計算
    Heff = H - (im/2) * sum(L'*L for L in Ls)
    
    # 時間発展の計算
    ts = 0:dt:tmax
    ψs = [copy(ψ0)]
    
    for t in ts[2:end]
        # ノルムの減少を計算
        norm_decrease = sum(abs2.(L*ψs[end]) for L in Ls)
        
        # ジャンプ確率の計算
        jump_prob = norm_decrease * dt
        
        if rand() < jump_prob
            # 量子ジャンプ
            ψ = normalize(sum(L*ψs[end] for L in Ls))
        else
            # 連続時間発展
            ψ = normalize(exp(-im*Heff*dt) * ψs[end])
        end
        push!(ψs, ψ)
    end
    
    return ts, ψs
end
```

## 使用例

### 1. ED法の使用例（example1_ED.jl）

```julia
include("../src/ED.jl")
include("../src/model.jl")

# パラメータの設定
params = Parameters(Δ=0.0, Γ=1/6, Ω=1.0)

# 初期状態の設定
ρ0 = sparse([0.0 0.0; 0.0 1.0])

# 時間発展の計算
ts = range(0, 40, length=1000)
H = makeHamiltonian(params)
L = makeLindbladian(params)
ρs = time_evol_ED(H, [1.], [L], ρ0, ts)
```

### 2. MCWF法の使用例（example2_MCWF.jl）

```julia
include("../src/MCWF.jl")
include("../src/model.jl")

# パラメータの設定
params = Parameters(Δ=0.0, Γ=1/6, Ω=1.0)

# 初期状態の設定
ψ0 = zeros(ComplexF64, 2)
ψ0[2] = 1.0

# 時間発展の計算
dt = 0.01
tmax = 20
ts, ψs = time_evol_trajectory(ψ0, dt, tmax, H, [L])
```

### 3. 両手法の比較（example3_comparison.jl）

```julia
include("../src/ED.jl")
include("../src/MCWF.jl")
include("../src/model.jl")

# パラメータの設定
params = Parameters(Δ=0.0, Γ=1/6, Ω=1.0)

# ED法による計算
ρs = time_evol_ED(H, [1.0], [L], ρ0, ts)
Pe_ED = [real(ρ[1,1]) for ρ in ρs]

# MCWF法による計算（1000サンプル）
n_samples = 1000
all_Pe = zeros(Float64, length(ts), n_samples)
for i in 1:n_samples
    _, ψs = time_evol_trajectory(ψ0, dt, tmax, H, [L])
    all_Pe[:,i] = [abs(ψs[j][1])^2 for j in 1:length(ts)]
end
```

## 結果の解釈

1. **ED法の結果**：
   - 正確な時間発展が得られる
   - 計算コストが高い（特に大規模系で）
   - 密度行列の完全な情報が得られる

2. **MCWF法の結果**：
   - 個々の量子軌道は確率的
   - アンサンブル平均で正確な結果に収束
   - 計算コストが低い（特に大規模系で）
   - 物理的な解釈が直感的

## まとめ

このチュートリアルでは、散逸のある量子系のシミュレーション手法として、ED法とMCWF法を実装しました。両手法にはそれぞれ長所・短所があり、問題に応じて適切な手法を選択することが重要です。

- ED法：小規模系、正確な結果が必要な場合
- MCWF法：大規模系、物理的な解釈が重要な場合

## 参考文献
1. H.-P. Breuer and F. Petruccione, "The Theory of Open Quantum Systems" (Oxford University Press, Oxford, 2002).
2. A. J. Daley, "Quantum trajectories and open many-body quantum systems", Advances in Physics 63, 77-149 (2014).

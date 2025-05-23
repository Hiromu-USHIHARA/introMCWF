このチュートリアルでは，開放量子系の数値シミュレーション手法である厳密対角化（exact diagonalization; ED）とモンテカルロ波動関数法（Monte Carlo wave function method; MCWF）を導入し，Juliaで実装します．
今回はNo. 1として，扱う問題の導入と厳密対角化による数値計算をおこないます．

## 扱う問題

このチュートリアルで扱うのは開放量子系と呼ばれるタイプの物理系です．
開放量子系は着目する系が外部自由度と相互作用しているものを指します．
例えば，熱浴と接触している物質や電磁場と相互作用する原子は開放量子系の例です．

開放量子系の記述として広く用いられているものに，Gorini--Kossakowski--Sudarshan--Lindblad方程式（GKSL方程式）があります^[適当な仮定のもとで開放量子系のダイナミクスがGKSL方程式に従うことが示されます．][1, 2, 3]:

$$
\frac{d\rho}{dt} = -i \left[ H, \rho(t) \right]_{-} + \sum_i \left( L_i \rho L_i^\dagger - \frac{1}{2} \left[ L_i^\dagger L_i, \rho \right]_{+} \right).
$$

ここで，$\rho$は着目する部分系の密度行列，$H$は着目系のHamiltonian，$L_i$はLindbladianと呼ばれる演算子で，$[A, B]_\pm=AB\pm BA$は（反）交換子です．

以下では，GKSL方程式で記述される開放量子系に注目して，その数値計算手法を解説します．
特に具体例として，以下で定義される電磁場と結合した二準位系を扱います[4]:

$$H=-\dfrac{\Omega}{2}\sigma_{\mathrm x}-\Delta\sigma_+\sigma_-,$$
$$L=\sqrt\Gamma\sigma_-.$$

## モデルの定義

まず，`model.jl`というファイルを作成して，２つの手法に共通のモデル定義を行います．

```julia
using LinearAlgebra, SparseArrays

mutable struct Parameters
    Δ::Float64
    Γ::Float64 
    Ω::Float64
    
    function Parameters(;Δ::Float64=0., Γ::Float64=1/6, Ω::Float64=1.0)
        new(Δ, Γ, Ω)
    end
end

σₓ = sparse(ComplexF64[0. 1.; 1. 0.])
σ₊ = sparse(ComplexF64[0. 0.; 1. 0.])
σ₋ = sparse(ComplexF64[0. 0.; 0. 1.])

function makeHamiltonian(p::Parameters)
    H = -p.Ω/2. * σₓ - p.Δ*σ₊*σ₋
    return H
end

function makeLindbladian(p::Parameters)
    L = √p.Γ * σ₋
    return L
end

function steadyPe(p::Parameters)
    Pe=abs(p.Ω)^2/4 / (p.Δ^2+p.Γ^2/4+abs(p.Ω)^2/2)
    return Pe
end

export Parameters, makeHamiltonian, makeLindbladian, steadyPe
```

## 数値解法１：厳密対角化（ED）

GKSL方程式は密度行列に関する方程式なので一見複雑に見えますが，数学的には1階の線型微分方程式にすぎないので（計算コストを無視すれば），係数行列の対角化によって解くことができます．
これが厳密対角化（ED）の基本的な方針です．

### Liouville space

注意すべき点として，GKSL方程式の係数"行列"は密度行列という演算子をその時間微分にマップするものです（超演算子と呼びます）．
したがって，対角化を実行するためには，超演算子を行列の形で書き換えてやるという準備が必要になります．
これを実現するのがLiouville spaceの方法です[5]．

Liouville spaceは密度行列の要素を1次元に並び替えて得られる線型空間のことです:

$$\rho=\begin{pmatrix}\rho_{11}&\rho_{12}&\cdots\\\rho_{21}&\rho_{22}\\\vdots&&\ddots\end{pmatrix}\longmapsto\left|\left.\rho\right>\right>=\begin{pmatrix}\rho_{11}\\\rho_{12}\\\vdots\\\rho_{21}\\\vdots\end{pmatrix}.$$

この変換に伴って超演算子の行列表現が得られます．
例えば，GKSL方程式に現れる$L_i\rho L_i^\dag$という項は，

$$L_i\rho L_i^\dag\longmapsto L_i\otimes L_i^\ast\left|\left.\rho\right>\right>$$

となります．
また，（反）交換子については，次のように書けることが知られています:

$$\left[A, B\right]_\pm\longmapsto [[A, 1]]_\pm\left|\left.B\right>\right>.$$

ここで超（反）交換子を

$$[[A, B]]_\pm:=A\otimes B^\top\pm B\otimes A^\top$$

によって定義しました．

以上の結果から，厳密対角化によってGKSL方程式を解くには，次の行列を数値対角化すれば良いということになります:

$$-i[[H, 1]]_- + \sum_i\left(L_i\otimes L_i^\ast - \dfrac12[[L^\dag L_i, 1]]_+\right).$$


### Juliaによる実装

さて，以上の方法をJuliaで実装していきましょう．

まず，Liouville spaceの方法の中心となる，$1$次元と$2$次元の間の配列の変換を実装していきます．


```julia
module sparse_transform_functions
    # 2Dのindex (m, n)に対し, M=N(m-1)+(n-1)+1を返す
    using SparseArrays
    function _2Dind_to_1Dind(m, n, N) # N=2^nsite
        return N*(m-1)+n
    end
    # 1Dのindex Mに対し, m, nを返す
    function _1Dind_to_2Dind(M, N)
        return (M-1)÷N + 1, (M-1)%N + 1
    end
    # %%
    # 2Dのmatrixに対し, 1Dの配列を返す
    function _2Dmat_to_1Darray(O_mat, N)
        result=spzeros(ComplexF64, N^2)
        for II ∈ 1:N^2
            m, n = _1Dind_to_2Dind(II, N)
            result[II] = O_mat[m, n]
        end
        return result
    end
    # 1Dの配列に対し, 2Dのmatrixを返す
    function _1Darray_to_2Dmat(O_array, N)
        result = spzeros(ComplexF64, N, N)
        for II ∈ 1:N^2
            m, n = _1Dind_to_2Dind(II, N)
            result[m, n] = O_array[II]
        end
        return result
    end
end
```

次に，Liouville spaceの表現を用いた行列表現及びその対角化による定常状態/ダイナミクス計算を実装します．

```julia
module sparse_Liouville_space
    using LinearAlgebra, KrylovKit, SparseArrays
    using Main.sparse_transform_functions
    ⊗(x, y)=kron(x, y)
    # %%
    function super_commutator(X, Y)
        X⊗transpose(Y) - Y⊗transpose(X)
    end
    #
    function super_anti_commutator(X, Y)
        X⊗transpose(Y) + Y⊗transpose(X)
    end
    #
    function dissipator_in_Liouville_space(γs, ℒs)
        @assert length(γs)==length(ℒs) "γs ($γs) & ℒs ($ℒs) must have same length"
        n = length(γs)
        𝕀 = sparse(I, size(ℒs[1]))
        #
        L = ℒs[1]
        result = γs[1]*(L⊗conj(L) - super_anti_commutator(L' * L, 𝕀)/2.)
        for i ∈ 2:n
            L = ℒs[i]
            result += γs[i]*(L⊗conj(L) - super_anti_commutator(L' * L, 𝕀)/2.)
        end
        return result
    end
    #
    function unitary_in_Liouville_space(H)
        𝕀 = sparse(I, size(H))
        return -im*super_commutator(H, 𝕀)
    end
    #
    function GKSL_ED_in_Liouville_space(H, γs, ℒs)
        𝔘 = unitary_in_Liouville_space(H)
        𝔇 = dissipator_in_Liouville_space(γs, ℒs)
        𝔏 = 𝔘 + 𝔇
        # es, vs = eigen(𝔏)
        es, vs, _ = eigsolve(𝔏, 3, :LR)
        @show es
        return es, vs
    end
    #
    function GKSL_NESS_in_Liouville_space(H, γs, ℒs)
        es, vs = GKSL_ED_in_Liouville_space(H, γs, ℒs)
        # @show inds_NESS = findall(x -> abs(x)<=1.0e-10, es)
        # result = transform_functions._1Darray_to_2Dmat(sum(vs[:, i] for i ∈ inds_NESS), size(H)[1]) # ρ_NESS
        result = sparse_transform_functions._1Darray_to_2Dmat(vs[1], size(H, 1)) # ρ_NESS
    end
    # 
    function time_evol_ED(H, γs, Ls, ρ0, ts)
        ℒ = unitary_in_Liouville_space(H) + dissipator_in_Liouville_space(γs, Ls)
        esL, vsL = eigen(Matrix(ℒ))
        #
        Vinv_ρ0 = vsL \ sparse_transform_functions._2Dmat_to_1Darray(ρ0, size(H, 1))
        @show sum(Vinv_ρ0[i]*vsL[:,i] for i ∈ 1:size(vsL, 2)) ≈ sparse_transform_functions._2Dmat_to_1Darray(ρ0, size(H, 1))
        ρs = [sparse_transform_functions._1Darray_to_2Dmat(sum(exp(esL[i]*t)* Vinv_ρ0[i] * vsL[:,i] for i ∈ 1:(size(H,1)^2)), size(H, 1)) for t ∈ ts]
        return ρs
    end
end
```

以下のコードを実行して，解析解と比較します．

```julia
include("../src/ED.jl")
include("../src/model.jl")

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
```

出力結果を下に示します．
振動しながら定常状態に緩和していくダイナミクスが観察できます．

![ED Results](assets/fig1_ED.png)


## No. 1のまとめ

今回は開放量子系を記述するGKSL方程式を導入し，厳密対角化によるダイナミクスの可視化を行いました．
次回の記事では，MCWF法による計算や２つの手法の比較を行う予定です．


## 参考文献
1. G. Lindblad, "On the generators of quantum dynamical semigroups", Commun. Math. Phys. 48, 119-130 (1976).
2. V. Gorini, A. Kossakowski, and E. C. G. Sudarshan, "Completely positive dynamical semigroups of N-level systems", J. Math. Phys. 17, 821 (1976).
3. H.-P. Breuer and F. Petruccione, "The Theory of Open Quantum Systems" (Oxford University Press, 2002).
4. A. J. Daley, "Quantum trajectories and open many-body quantum systems", Adv. in Phys. 63, 77-149 (2014).
5. J. A. Gyamfi, "Fundamentals of quantum mechanics in Liouville space", Eur. J. Phys. 41, 063002 (2020).

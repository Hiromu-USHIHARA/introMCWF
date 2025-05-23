この記事では，[前回](https://zenn.dev/hiromu_ushihara/articles/db67c435b9b95b)に続いて，開放量子系の数値計算手法について扱います．
特に今回はMonte Carlo Wave Function (MCWF)法の導入・実装と厳密対角化（ED）との比較を行います．

> このチュートリアルは[Zennでも公開](https://zenn.dev/hiromu_ushihara/articles/c6ef07f16666ee)していますので，数式がうまく表示されない場合はそちらをご覧ください．

## 扱う問題（再確認）

このチュートリアルでは，Gorini--Kosakkowski--Sudarshan--Lindblad (GKSL)方程式[1, 2, 3]

$$
\frac{d\rho}{dt} = -i \left[ H, \rho(t) \right]_{-} + \sum_i \left( L_i \rho L_i^\dagger - \frac{1}{2} \left[ L_i^\dagger L_i, \rho \right]_{+} \right).
$$

で記述される開放量子系を扱います．
特に具体例として，以下で定義される電磁場と結合した二準位系を扱います[4]:

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

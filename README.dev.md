# ローカル開発とRSM・PCM検証

開発ブランチは `dev/rsm-design-diagnostics` です。初期の重点は基本RSMにおける疎・不均衡な
採点計画の識別性と推定精度です。再現実験は隣接する `../mfrm-replication-bayes` から実行します。

## 確認手順

```sh
Rscript --vanilla tools/test-local.R
Rscript --vanilla ../mfrm-replication-bayes/scripts/run-small-study.R
```

通常のテストはMCMCを起動しません。既存の全モデル統合テストを明示的に実行する場合は
`BAYESMFRM_RUN_STAN=true` を設定します。今回の小規模実験は基本RSMの推定と出力を確認します。
`renv.lock` は現在の外部依存環境を記録し、必要なときに `renv::restore()` で復元します。
CmdStanのインストールは別途必要です。小規模実験では共有PCHの書換えを避けるため
`cpp_options = list(PRECOMPILED_HEADERS = "false")` を指定します。

## 採点計画とモデル

```r
library(BayesMFRM)
plan <- expand.grid(person = 1:20, item = 1:4, rater = 1:4)
design <- check_design_bmfrm(score ~ person + item + rater, plan, K = 6)
print(design)
prior <- prior_predict_bmfrm(score ~ person + item + rater, plan, K = 6,
                             ndraws = 200, seed = 42)
sim <- simulate_bmfrm(score ~ person + item + rater, plan, K = 6, seed = 43)
fit <- bmfrm(score ~ person + item + rater, sim$data, K = 6,
              design_check = "error", seed = 44)
bayesplot::pp_check(fit, group = "rater", type = "all")
residual_cor_criteria(residuals(fit), unit_facets = c("person", "rater"))
```

隣接カテゴリの対数オッズは `theta[person] - item[item] - rater[rater] - tau[step]` です。
受験者能力は中心化せず、その他の主効果と共通閾値に和ゼロ制約を適用します。
計画検査はこの制約を適用した予測子行列に直交基底を使い、SVDによる階数と零空間を返します。
共通閾値の自由度 `K - 2` は予測子の自由度とは独立です。連結性だけでは識別性を保証しません。

検査はRSM・PCMと2因子交互作用を対象にします。任意のStanファイルの識別性は認証しません。
代数的識別性、推定精度、極端評定による最尤推定の発散は別の問題です。
`observed_only = TRUE` で欠測評定を除いた計画を検査できます。明示的な因子水準は未観測でも保持します。
初期実装は密行列を使い、過大な計画は計算前にエラーにします。

事前生成は `normal()`・`student_t()`・`double_exponential()` の数値引数をサポートし、Stanと同じ事前分布・制約を使います。
中心化後の事前分散は元の独立事前の分散とは異なります。事前予測は妥当性を自動判定するものではありません。

残差相関は記述統計です。複数評価者がいる場合の対応単位を `unit_facets` で明示します。
受験者×項目単位に平均する場合は `aggregate = "mean"` を明示してください。
公平得点の区間は同じ事後サンプル内で集約し、観測間の共分散を保持します。
ただし現在の公平得点は各受験者が評価された項目構成での平均で、共通項目への標準化は未実装です。

キャッシュの署名にはモデル・データ・サンプリング設定・CmdStan版を含めます。
署名がない、設定が異なる、チェーンのCSVが失われたキャッシュは再利用しません。
開発時の検証結果は隣接リポジトリの `reports/validation.md` に記録しています。

## 開発段階と今後の課題

### 引数によるモデル指定（開発版 0.1.0.9001）

公開コンストラクタ `rsm_model()` と `model` 引数を廃止し、既存の `family`・`priors` と
ファセット名を受け取る `hierarchical` 引数に整理しました。新しい公開関数は追加しません。
最小の指定は `bmfrm(..., family = "rsm", hierarchical = "person")` です。
既定の `hierarchical = character(0)` は従来の独立事前分布を維持します。

```r
priors <- c(
  prior("normal(0, 2.5)", class = "mean", facet = "person"),
  prior("lognormal(0, 0.5)", class = "sd", facet = "person"),
  prior("normal(0, 1)", class = "item"),
  prior("normal(0, 1)", class = "rater"),
  prior("normal(0, 1)", class = "tau")
)
sim <- simulate_bmfrm(score ~ person + item + rater, plan, K = 6,
                     hierarchical = "person", priors = priors, seed = 42)
fit <- bmfrm(score ~ person + item + rater, sim$data, K = 6,
             hierarchical = "person", priors = priors, design_check = "error", seed = 43)
```

能力は `theta = mu_theta + sigma_theta * z`、`z ~ normal(0, 1)` です。
`prior(0, "mean", facet = "person")` や `prior(2, "sd", facet = "person")` ならその母集団
パラメータを固定します。分布・固定値とも同じ `priors` で指定し、生成・推定が共通の内部検証を使います。
独立能力の `theta` 事前と階層化された受験者の `mean` / `sd` を混在させるとエラーになります。
明示的な `priors` は必要な対象をすべて指定してください。
既存の `get_default_priors(..., hierarchical = "person")` でも完全な一覧を取得できます。

**事前SDの意味は階層化の有無で変わりません。** `prior()` の既定値 `scale = "raw"` は
中心化前のSDを表します。J水準の正規効果を和ゼロにすると、各水準の分散は `sd^2 * (1 - 1/J)` です。
`scale = "marginal"` を明示した `normal(0, sd)` では、各水準の周辺分散が `sd^2` になります。
いずれも1水準の効果と2値RSMのステップは0です。
Gaussianブロックでは前者を `sd * H * z`、後者を `sd * sqrt(J/(J-1)) * H * z` と実装します。
基底Hの選択は内部処理です。旧来の独立・raw指定の経路はStudent-t事前を含めて維持します。
階層化またはmarginal指定の新経路は、Gaussian RSM・PCMと2因子交互作用に対応します。

旧 `rsm_model()` と同じ事前分布に移行するには、母集団のmean/sdを指定したうえで、
item・rater・tauそれぞれの `prior()` に `scale = "marginal"` を明示します。
`hierarchical = "person"` だけでは旧コンストラクタの周辺SD設定を再現しません。
変更前コミット `b0d059c` から保存した生成値・Stan入力・生成コードと照合する回帰テストがあります。

ファセットが10個に増えても、式に列名を列挙し、事前分布を既存の `facet_<name>` で指定します。
ファセットごとの専用引数は増やしません。母集団分布の指定先は `class` と `facet` の組で識別します。
主効果の階層化対応先は受験者 `person` のみです。`hierarchical = c("person", "rater")` は
未実装としてエラーにします。評価者などの和ゼロ制約と母集団パラメータの識別性を検証してから拡張します。

MCMC診断は自由座標と復元したパラメータを確認し、構造的に固定された量を除外します。
計画の階数検査は直接の能力・facet効果の条件付き予測子を対象とし、母集団平均・SDの
実用的な識別や精度を保証しません。明示された未観測因子水準は全工程で保持します。

```sh
# 独立R数式とStanの対数密度差・解析勾配を固定点で比較（MCMCなし）
Rscript --vanilla tools/check-rsm-stan.R
```

公開APIは既存の入口と役割が重ならない最小限の引数で拡張します。
内部の仕様共有のためだけに公開コンストラクタを追加しません。
既存設定の意味、生成・推定の一致、固定点パリティを検証してから実推定を確認します。
小規模な動作確認から一般的な校正や推定精度改善は主張しません。

### Step priorとPCMの識別制約（開発版 0.1.0.9002）

Stepにも事前分布の設計が必要です。既存の `prior(..., class = "tau")` を使い、
公開関数・事前クラスは増やしていません。RSMでは全水準で共有するstep、PCMでは
`threshold_facet` の各水準に固有のstepを指定します。

```r
ps <- c(
  prior("normal(0, 2)", "theta"),
  prior("normal(0, 1)", "item"),
  prior("normal(0, 1)", "rater"),
  prior("normal(0, 0.7)", "tau", scale = "marginal")
)
sim <- simulate_bmfrm(score ~ person + item + rater, plan, K = 4,
  family = "pcm", threshold_facet = "item", priors = ps, seed = 45)
pp <- prior_predict_bmfrm(score ~ person + item + rater, plan, K = 4,
  family = "pcm", threshold_facet = "item", priors = ps, ndraws = 200, seed = 46)
fit <- bmfrm(score ~ person + item + rater, sim$data, K = 4,
  family = "pcm", threshold_facet = "item", priors = ps, design_check = "error", seed = 47)
```

PCMでは隣接対数オッズを `theta[p] - item[i] - rater[r] - tau[i,k]` とし、
**各項目内で `sum(tau[i,]) = 0`** とします。`item[i]` は全体難易度、`tau[i,k]` は
その難易度からのstep偏差です。絶対的な隣接閾値を調べる場合は、同じ事後drawの
`item[i] + tau[i,k]` を合成します。要約後の平均・区間を足してはいけません。
`summary()` は行列添字と元のファセット水準名を表示します。

`S = K-1 > 1` とすると、raw Gaussian priorの各stepの周辺分散は `sd^2*(1-1/S)`、
marginal指定では `sd^2`、異なるstep間の共分散はそれぞれ `-sd^2/S` と `-sd^2/(S-1)` です。
この尺度は項目数に依存しません。項目間のstepブロックは独立です。
2値PCMではstep自由度が0になり、1項目の多値PCMでは `K-2` 自由度が残ります。
raw Student-t priorも項目内で中心化しますが、Gaussian用の周辺SD補正は適用しません。
和ゼロ化されるブロックでは、全要素に共通の事前位置を変えても中心化で相殺されます。

このpriorはstep偏差を0へ正則化します。項目間で共通のstepパターンへ部分プーリングする
階層prior、step尺度のハイパーprior、順序制約は未実装です。順序制約はカテゴリの順序と
異なる追加仮定であり、デフォルトでは設けません。尺度の妥当性は、対象計画の事前予測で
カテゴリの出現割合や極端評定を確認し、複数の尺度で感度分析してください。

**旧PCMから推定対象が変わります。** 旧実装はstepを項目間で中心化していたため、
和ゼロの任意の `c[i]` について `item[i] + c[i]` と `tau[i,k] - c[i]` が同じ尤度を与えました。
さらに、共通stepパターンも制限していました。proper priorで事後分布を定義できることと、
尤度が識別可能なことは別です。新旧PCMの結果を同じモデルの再推定として混在させないでください。
生成Stanコードが変わるため、旧fitキャッシュは再利用されません。RSMの既存設定は維持します。

PCMの計画検査は、主効果の階数に観測されたthreshold水準ごとの `K-2` を加えます。
未観測水準のstep自由度は `step_nullity`・`unobserved_step_levels` に報告します。
`rank`・`null_space`・`condition_number` は主効果予測子の値であり、全隣接対数オッズの
階数は `adjacent_logit_rank` です。カテゴリの観測頻度による精度・実用的識別性は別途確認します。

`threshold_facet = "rater"` など、調整対象にstepもある場合の公平得点は、参照stepを
定義するまで `fair_scores()` で明示的に停止します。項目別stepの場合は従来どおり、
観測された項目構成での得点を返します。異なる採点計画間の標準化はまだ行いません。

### 交互作用の階層priorと診断（開発版 0.1.0.9003）

2因子交互作用の指定には既存の式・`hierarchical`・`priors` を使います。
公開関数やpriorクラスは追加していません。交互作用名の成分順は照合時に同一視し、
生成Stanの行列次元はデータ側の水準数から宣言します。

```r
f <- score ~ person + item + rater + item:rater
h <- c("person", "item:rater")
ps <- c(
  prior("normal(0, 1)", "mean", facet = "person"),
  prior("lognormal(0, 0.35)", "sd", facet = "person"),
  prior("normal(0, 0.7)", "item", scale = "marginal"),
  prior("normal(0, 0.7)", "rater", scale = "marginal"),
  prior("normal(0, 0.7)", "tau", scale = "marginal"),
  prior("normal(0, 0.5)", "sd", facet = "rater:item", scale = "marginal")
)
check_design_bmfrm(f, plan, K = 4)
sim <- simulate_bmfrm(f, plan, K = 4, hierarchical = h, priors = ps, seed = 51)
pp <- prior_predict_bmfrm(f, plan, K = 4, hierarchical = h, priors = ps, seed = 52)
fit <- bmfrm(f, sim$data, K = 4, hierarchical = h, priors = ps,
             design_check = "error", seed = 53)
summarise_bias(fit, "item:rater", delta = 0.3)
```

`class = "sd"` は、階層化した交互作用全体の尺度に対するpriorです。
`normal(0,0.5)` は正の範囲に制限され、half-normalになります。この0.5はhalf-normalの
元の正規分布の尺度であり、交互作用のSDを0.5に固定する指定ではありません。
`student_t(df,0,s)`・`cauchy(0,s)` の正の範囲、`lognormal(m,s)` も指定できます。
`prior(0.5,"sd",facet="item:rater")` は尺度の固定値です。
一方、階層化しない交互作用は従来どおり `class = "bias"` で効果の分布を指定します。
同じ交互作用に両方を指定すること、和ゼロ効果に母集団平均を追加することはエラーになります。

交互作用行列Bの各行・各列を、**宣言された水準について等重みで和ゼロ**にします。
正規経路は `B = sigma * H1 * Z * t(H2)`、Zの各要素が標準正規です。
自由度は `(J1-1)*(J2-1)`。raw指定ならHiは正規直交基底で、各セルの条件付き分散は
`sigma^2*(1-1/J1)*(1-1/J2)` です。marginal指定では両基底を補正し、条件付き分散を
`sigma^2` に揃えます。いずれもセル間は独立ではありません。どちらかが1水準ならBは0で、
尺度を推定する指定があっても、その尺度への情報はpriorのみです。
基底の選択は内部処理であり、等方的な正規priorは直交基底の取り替えに依存しません。

正負の符号は予測子に `+B` として入ります。正の交互作用は、他の効果が同じなら高い評定を
取りやすくする、モデル上の相対的な偏りです。
固定値で生成する場合は、式が定める成分順の行列（例 `itemrater`）を渡します。
任意の2因子名や複数の交互作用に対応します。名称がStan変数として衝突する場合と3因子以上は停止します。

rawの非正規priorは各生セルに置いて二重中心化します。中心化後の各セルが同じStudent-tや
ラプラス分布になるとは主張しません。階層化・marginal経路では、主効果・step・条件付き交互作用は
Gaussianに限定します。LKJや馬蹄prior、主効果の評価者階層化は未実装です。

`check_design_bmfrm()` は交互作用の対比も含めた階数と零空間を返します。
未観測セルも `interaction_levels` に残り、観測数と制約写像の行番号を持ちます。
グラフの連結性だけでは交互作用を含むモデルの識別性は保証されません。
構造的識別性の検査は効果を直接の座標として扱い、尺度の実用的な推定精度を保証しません。

`summarise_bias()` は生成されたパラメータ名・成分順・元のラベルを使って全セルを要約します。
`delta` は隣接対数オッズでの実質的な差です。方向確率と `P(|B| > delta | data)` を返し、
未観測数・尤度で識別できるかも併記します。識別できないセルは事後要約を保持したまま
`not_likelihood_identified` とし、バイアス検出に数えません。95%区間と確率閾値によるフラグは
多重比較の誤検出率制御ではありません。旧 `bias[I,R]` という生CmdStan出力には
既存の `summarise_bias_ir()` を使います。

公平得点ではrater/interlocutorの主効果と、それらを含む交互作用を除きます。
それ以外の交互作用は保持します。公平得点は依然として観測された項目構成に依存します。

```sh
Rscript --vanilla tools/check-bias-stan.R
```

### 続く拡充

- 共通stepと水準別偏差の階層化を、独立PCM・RSMとの関係と識別制約を定義して検討する。
- 交互作用の検出対象と、セル・評価者・項のどこに疎性を仮定するかを定義して、構造化priorを比較する。
- 採点計画と事前分布ごとのSBC、事前感度、残差構造の事後予測チェックを追加する。
- 共通項目に標準化した公平得点と、大規模計画用の疎行列による階数検査を追加する。
- R CMD checkの静的解析NOTEを解消し、全モデルの統合テストを拡充する。

PRはモデル定義、実装、回帰テスト、利用者向け説明を同じ単位で更新します。

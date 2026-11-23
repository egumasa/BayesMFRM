# ローカル開発とRSM検証

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

検査は加法的RSMに限定します。PCM・交互作用・任意のStanファイルの識別性は認証しません。
代数的識別性、推定精度、極端評定による最尤推定の発散は別の問題です。
`observed_only = TRUE` で欠測評定を除いた計画を検査できます。明示的な因子水準は未観測でも保持します。
初期実装は密行列を使い、過大な計画は計算前にエラーにします。

事前生成は `normal()` と `student_t()` の数値引数をサポートし、Stanと同じ生パラメータの生成・中心化を行います。
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
階層化またはmarginal指定の新経路は、現在は加法的Gaussian RSMに限定します。

旧 `rsm_model()` と同じ事前分布に移行するには、母集団のmean/sdを指定したうえで、
item・rater・tauそれぞれの `prior()` に `scale = "marginal"` を明示します。
`hierarchical = "person"` だけでは旧コンストラクタの周辺SD設定を再現しません。
変更前コミット `b0d059c` から保存した生成値・Stan入力・生成コードと照合する回帰テストがあります。

ファセットが10個に増えても、式に列名を列挙し、事前分布を既存の `facet_<name>` で指定します。
ファセットごとの専用引数は増やしません。母集団分布の指定先は `class` と `facet` の組で識別します。
現在の階層化対応先は受験者 `person` のみです。`hierarchical = c("person", "rater")` は
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

### 続く拡充

- PCMの項目効果と項目別閾値の分離を修正し、独立した識別性テストを追加する。
- 交互作用の生成・要約で使うパラメータ名と出力の契約を統一する。
- 採点計画と事前分布ごとのSBC、事前感度、残差構造の事後予測チェックを追加する。
- 共通項目に標準化した公平得点と、大規模計画用の疎行列による階数検査を追加する。
- R CMD checkの静的解析NOTEを解消し、全モデルの統合テストを拡充する。

PRはモデル定義、実装、回帰テスト、利用者向け説明を同じ単位で更新します。

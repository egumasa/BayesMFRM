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

### 第2段階: 明示的な等距RSM仕様

```r
model <- rsm_model(
  ability_mean = "normal(0, 2.5)",
  ability_sd = "lognormal(0, 0.5)",
  facet_sd = c(item = 1, rater = 1), step_sd = 1
)
sim <- simulate_bmfrm(score ~ person + item + rater, plan, K = 6,
                     model = model, seed = 42)
fit <- bmfrm(score ~ person + item + rater, sim$data, K = 6,
             model = model, design_check = "error", seed = 43)
```

`model = NULL` は従来の生パラメータ事前分布を維持します。新しい仕様は加法的RSMのみで、
`priors` 引数と併用しません。`ability_mean = 0, ability_sd = 2` のような数値を指定すれば
対応する母集団パラメータを固定できます。文字列で指定すると平均は正規事前、SDは対数正規事前で推定します。

新仕様では `theta = mu_theta + sigma_theta * z`、`z ~ normal(0, 1)` です。
J水準の効果は `x = sd * sqrt(J / (J - 1)) * H * z` とし、Hは正規直交和ゼロ基底です。
したがって全水準の周辺事前分散が `sd^2` になります。1水準の効果と2値RSMのステップは恒等的に0です。
これは従来の `raw - mean(raw)` と同じ数値SDを代入して同じ事前分布になる変更ではありません。
生成と推定のモデル指定を揃え、ラベル順序・共分散・固定点対数密度・勾配で検証します。

新仕様のMCMC診断は自由座標と復元したモデルパラメータを確認し、構造的に固定された量を除外します。
計画の階数検査は直接の能力・facet効果の条件付き予測子を対象とし、母集団平均・SDの
実用的な識別や精度を保証しません。明示された未観測因子水準も生成・推定・計画検査で保持します。

```sh
# 独立R数式とStanの対数密度差・解析勾配を固定点で比較（MCMCなし）
Rscript --vanilla tools/check-rsm-stan.R
```

この段階の終了条件は、(1) 既存テストを保持、(2) ラベル置換と事前共分散を確認、
(3) 母集団の推定／固定・2値／多値・1水準の固定点パリティ、(4) 全試行を保持した小規模実推定です。
一般的な校正や精度改善を確認したという主張には、別途十分なSBC・回復実験が必要です。

### 続く拡充

- PCMの項目効果と項目別閾値の分離を修正し、独立した識別性テストを追加する。
- 交互作用の生成・要約で使うパラメータ名と出力の契約を統一する。
- 採点計画と事前分布ごとのSBC、事前感度、残差構造の事後予測チェックを追加する。
- 共通項目に標準化した公平得点と、大規模計画用の疎行列による階数検査を追加する。
- R CMD checkの静的解析NOTEを解消し、全モデルの統合テストを拡充する。

PRはモデル定義、実装、回帰テスト、利用者向け説明を同じ単位で更新します。

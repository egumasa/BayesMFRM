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

## 今後の課題

- PCMの項目効果と項目別閾値の分離を修正し、独立した識別性テストを追加する。
- 交互作用の生成・要約で使うパラメータ名と出力の契約を統一する。
- 採点計画と事前分布ごとのSBC、事前感度、残差構造の事後予測チェックを追加する。
- 共通項目に標準化した公平得点と、大規模計画用の疎行列による階数検査を追加する。
- R CMD checkの静的解析NOTEを解消し、全モデルの統合テストを拡充する。

PRはモデル定義、実装、回帰テスト、利用者向け説明を同じ単位で更新します。

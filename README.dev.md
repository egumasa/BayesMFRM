# BayesMFRM 開発・検証ノート（v0.1.0.9007 現行仕様）

共同編集ブランチは `fix/audit-findings` です。Git管理されたチェックアウトを編集・検証・インストールの正本としてください。この作業領域では `local-work/BayesMFRM` が該当し、その親にあるコピーは独立した開発場所として扱いません。再現実験の別リポジトリは、このチェックアウトに隣接する `../mfrm-replication-bayes` に配置しています。過去の設計変更・破壊的変更の履歴は [NEWS.md](NEWS.md)、利用者向け入門ガイドは [README.md](README.md) および `vignette("getting-started", package = "BayesMFRM")`（[vignettes/getting-started.Rmd](vignettes/getting-started.Rmd)）を参照してください。

## テストと検証コマンド

```sh
# 1. 単体テスト（MCMCは起動しない）
Rscript --vanilla tools/test-local.R

# 2. Stan 対数密度・解析勾配の固定点パリティ検証（独立R実装と照合）
Rscript --vanilla tools/check-rsm-stan.R
Rscript --vanilla tools/check-bias-stan.R
Rscript --vanilla tools/check-extensions-stan.R

# 3. 全モデル統合テスト（実MCMCを実行する場合のみ明示指定）
BAYESMFRM_RUN_STAN=true Rscript --vanilla tools/test-local.R
```

共有プリコンパイル済みヘッダの書換えを避けるため、再現実験では `cpp_options = list(PRECOMPILED_HEADERS = "false")` を指定します。このチェックアウトの `renv.lock` に外部依存を固定しています。復元はチェックアウトのルートで `renv::restore()` を実行してください。CmdStan 本体は別途インストールが必要です。`Rscript --vanilla` は自動的なrenv有効化を省くため、上記の検証では実際に使われたライブラリと `sessionInfo()` も確認してください。

例示的なBayes対ML実験は `Rscript --vanilla simulations/bayesian_benefits.R` で実行できます。出力は `simulations/results/` に保存され、Gitには含めません。`docs/reproduce-ratings1.R` とこの実験スクリプトは、ルートから実行すると現在のソースを `pkgload` で読み込みます。

**0.1.0.9007の監査修正の検証（2026-11-27、macOS arm64 / R 4.6.1 / CmdStan 2.39.0）**

- 通常テスト: 1,107 assertions成功、失敗・エラー・警告0件。既存の実MCMC統合テスト3件は既定どおりスキップ。
- ビネットを含むソースビルドと `R CMD check --no-manual --no-tests`: `Status: OK`。単体テストは上記で別実行し、重複実行を省いた。
- 一時ライブラリへインストールしたパッケージで `example_usage.R` の2モデルを実MCMC推定し、全後処理を実行。fitの保存・再読込後も残差が一致し、保存Stan対R再計算のLOOの観測別ELPDの最大差は `4.9e-14`。
- 同梱カスタムStanの実推定で、使用していないR側の事前分布がメタデータに記録されないことを確認。推定済みキャッシュの再読込も確認。
- 再現実験スクリプトは初期化・現行テンプレートからのコード生成・Stan検査まで確認。保存済み実験結果の再計算や、大規模シミュレーション・他OSでの検証は行っていない。使用例と短いカスタムStan推定には収束警告があり、これらは実行確認であって推定精度の保証ではない。

## 現行モデル仕様と制約

- **線形予測子と和ゼロ制約**: 隣接カテゴリ対数オッズは RSM で `theta[person] - sum(facet) + sum(bias) - tau[k]`、PCM で `tau[threshold_facet, k]` とします。受験者能力 `theta` は中心化せず、その他の主効果・ステップ偏差・2因子交互作用（行・列）は宣言された水準について等重みで和ゼロ制約を課します。PCM のステップ偏差は **`threshold_facet` の各水準内で和ゼロ** です。
- **事前分布と尺度 (`scale`)**:
  - `scale = "raw"`（既定）：中心化前のSD。$J$ 水準の正規ブロックでは制約後の周辺分散が $\sigma^2(1 - 1/J)$、2因子交互作用では $\sigma^2(1 - 1/J_1)(1 - 1/J_2)$ になります。
  - `scale = "marginal"`：制約後の各水準・各セルの周辺分散を $\sigma^2$ に揃えます（Gaussian ブロックおよび単位分散正規化 Student-t 局所収縮 `inv_gamma(1.5, 0.5)`）。
  - `local_shrinkage = "horseshoe"`：分布固有のグローバル収縮尺度のため `scale = "raw"` のみ許可します。潜在セル $U$ に局所収縮をかけた後に二重中心化 $B = P_1 U P_2$ を行うため、セル間波及が生じます。
  - 非正規事前（`student_t`, `double_exponential`）と階層・`marginal` 指定の混在に対応します。
- **階層化 (`hierarchical`)**: `"person"`（`mean`/`sd`）、任意の和ゼロ主効果（`"rater"`, `"item"` 等、`class = "sd"`）、2因子交互作用（`class = "sd"` + 任意 `local_shrinkage`）、階層 PCM（`"tau"`：`tau_common` への部分プーリングと水準別偏差）に対応します。
- **識別子・`K` 検証**: 予約済み Stan 識別子（`theta`, `tau`, `*_raw`, `model` 等）の検査は `parse_bmfrm_formula()` および `resolve_bmfrm_spec()` の分岐前で統一実行します。`K` 省略時は `max(score)` から推定する旨のメッセージを出します。

## 診断・後処理・スケーラビリティ設計

- **計画識別性検査 (`check_design_bmfrm()`)**: 密行列上限を超える大規模計画（例: 1,500人・6,000観測以上）では `method = "auto"` が厳密な Schur 補行列（within-person TSQR＋Lanczos による全設計行列の exact 条件数）ソルバへ自動切替します。サイズ上限検査は基底確保より前に実行し、Schur モードでは Kronecker 基底を実体化せず `kron_rows()` / `apply_kron_null()` で計算します（片側1水準の自由度0交互作用では基底確保をスキップ）。
- **軽量保存とストリーミング後処理 (`save_gq = "minimal"`)**: 観測長 `N` の `generated quantities` を省略し、`R/posterior_quantities.R` の 2 パス softmax＋チャンク累積で `residuals()`, `facet_fit()`, `fair_scores()`, `bayesplot::pp_check()` を $O(D \cdot J)$ メモリで計算します。
- **内部評価契約 (`resolve_eval_contract()`)**: `spec$eval_contract` が生成モデルおよびバイト一致で検証済みの同梱ハードコードモデルのみに R 側再計算を許可します。未検証のカスタム Stan ファイルでは変数名が一致していても R 再計算を拒否し、Stan 保存済み生成量のみを参照します。
- **公平得点 (`fair_scores()`)**: `adjust_facets`、`standardize = TRUE`（巨大グリッド展開なしのチャンク累積、`max_grid_cells` / `max_output_rows` ガード、`person` 交互作用・受験者別 PCM 対応）、`likelihood_identified` および `identification_basis`（`conditional_likelihood` / `hierarchical_step_prior` / `unidentified_predictor_or_step`；自由ステップのない 2 値 PCM `K = 2` は事前依存扱いしない）を返します。
- **適合度・バイアス・キャッシュ**:
  - `facet_fit()` は各事後 draw 内で 4 次中心モーメントから $q_d$ と Wright & Masters (1982) Wilson-Hilferty ZSTD を算出します。`facet_infit_outfit()` / `add_fit_cri_bayes()` は裸の列名と文字列の双方を受理します。
  - `summarise_bias()` と `summarise_bias_ir()` は共通の `classify_bias_flags()` で単一の方向確率閾値 `prob_thresh` を適用し、`summarise_bias()` は同順位不変な `bfdr` を返します。
  - `fit_cmdstan_cached()` はモデル・データ・サンプリング設定・CmdStan版を署名に含め、追加 `...` なしの最小呼び出しにも対応します。
  - `print.bmfrm_fit()` / `summary.bmfrm_fit()` は計画識別性と MCMC 収束を分離して表示します。

## 今後の課題

- R CMD check の静的解析 NOTE の解消および実 MCMC を含む統合テストの拡充。
- 事前登録型の小規模・大規模 SBC および比較実験（`../mfrm-replication-bayes/scripts/`）の実実行と結果アーカイブ。
- 公開ドキュメントサイト（`pkgdown` 等）の整備。

# 一次資料の所在

このプロジェクトの方針は「パラメータや定数を記憶や推測で書かず、
必ず出典を書く」である。到達できているもの・いないものをここに一覧する。

## 到達できているもの

| # | 資料 | 所在 | 用途 |
|:--|:---|:---|:---|
| **S** | **Falcon 仕様書** *Falcon: Fast-Fourier Lattice-based Compact Signatures over NTRU*, **Specification v1.2 — 01/10/2020**、67 頁 | 本文はリポジトリに含めない（下記） | `[Spec]` タグの引用元。**パラメータの一次資料** |
| A | **Falcon 参照実装（C）** `Falcon-impl-20211101.zip` | `scripts/cref/`（vendor 済、MIT）。経緯は `scripts/cref/PROVENANCE.md` | `[C-ref]` タグの引用元。速度対決の相手 |
| B | **Python 参照実装** tprest/falcon.py | `scripts/pyref/`（vendor 済、MIT） | `[Py-ref]` タグの引用元。golden vector 生成 |
| C | **Pornin, Prest: More Efficient Algorithms for the NTRU Key Generation using the Field Norm** (IACR ePrint **2019/015**, 31 頁) | 本文はリポジトリに含めない（下記） | 降下アルゴリズムの一次資料 |

### S（仕様書）のどこに何があるか

`test/test_spec.jl` が該当箇所を**転記して検査している**ので、
この表は入口であって、実物はそちらにある。

| 箇所 | 内容 | 対応 |
|:---|:---|:---|
| **Table 3.3**（p.51） | パラメータ表。σ, σmin, σmax, ⌊β²⌋, 公開鍵長, 署名長 | `FALCON_512` / `FALCON_1024`。`spec_ref` はここを指す |
| (2.10) §2.6 | `q = 12·1024+1 = 12289`、「`k·2n+1` 型の最小の素数」 | `Q` |
| (2.11) §2.6 | `‖B‖_GS ≤ 1.17·√q` | `gram_schmidt_quality()` |
| (2.12) §2.6 | `σ_{f,g} = 1.17·√(q/2n)` | `SIGMA_FG_BASE` と `sigma_fg` |
| (2.13) §2.6 | `σ = (1/π)·√(log(4n(1+1/ε))/2)·1.17·√q`、`ε ≤ 1/√(Q_s·λ)` | `sigma`。**σmin とは別の ε**（#046） |
| (2.14) §2.6 | `β = τ·σ·√(2n)`, `τ = 1.1`、判定は `⌊β²⌋` | `sig_bound` |
| **Table 3.1**（p.41） | BaseSampler の pdt/cdt/RCDT（2^72 倍） | `RCDT`, `RCDT_PREC` |
| **Table 3.2**（p.44-45） | SamplerZ のテストベクタ 16 本 | `test_spec.jl`。**#022 の replay 規約はここで決着**（`reversed_chunks = true`） |
| Algorithm 1-18 | `splitfft`/`mergefft`/`HashToPoint`/`Keygen`/`NTRUGen`/`NTRUSolve`/`Reduce`/`LDL*`/`ffLDL*`/`Sign`/`ffSampling`/`BaseSampler`/`BerExp`/`SamplerZ`/`Verify`/`Compress`/`Decompress` | 各モジュール |
| §4 | 実装上の注記（浮動小数点、FFT/NTT、LDL 木、鍵生成、性能） | `docs/math/` |

**σmin だけは Table 3.3 に値しかなく、導出が書かれていない。**
`smoothing_eta` / `falcon_eps` は依然 `[derived]` である（#046）。

**秘密鍵のバイト長も Table 3.3 に無い**
（"Private key size (not listed above) is about three times that of a signature"）。
`privkey_bytes` は `[C-ref]` のままである。

### S（仕様書）と C（論文）をリポジトリに入れていない理由

MIT ではないので、全文を vendor しない。
**節番号で引用し、必要な箇所だけ短く引く**という形にしてある
（この文書と `docs/debug_log.md`、`docs/math/06_ntrugen.md` がそれ）。
入手先: 仕様書は <https://falcon-sign.info/>（提出パッケージ同梱の
`falcon.pdf`）、論文は <https://eprint.iacr.org/2019/015>。

### 論文のどこに何が書いてあるか（この実装に関係する範囲）

| 節 | 内容 | この実装での対応 |
|:---|:---|:---|
| §2.8, Algorithm 1 | `Reduce` の理想版。`k ← ⌊(Ff*+Gg*)/(ff*+gg*)⌉` を `k ≠ 0` の間くり返す | `main` ブランチの `babai_reduce` |
| §2.8 の地の文 | 「浮動小数点なので無限ループしうるが、**(F,G) のノルムが減らなくなった時点で抜ければよい**」 | 停止条件は仕様の `k = 0` に**限定されていない**。分岐版の根拠 |
| §2.9, Algorithm 2 | `ResultantSolver`（旧来の終結式による解法、O(n³)） | 実装していない（比較対象として言及のみ） |
| §3, Algorithm 3 | `TowerResultant`、体ノルムによる終結式 | `field_norm` |
| §4.2, Algorithm 4 | `TowerSolverR`（再帰版） | `ntru_solve` はこちらの形 |
| §4.3, Algorithm 5 | `TowerSolverI`（反復版、省メモリ） | C 参照実装はこちら |
| §5.2 | RNS / CRT / NTT、31 bit ワード | `docs/math/06_ntrugen.md` 6.9b |
| §5.3 | 二進 GCD を 31 段まとめて進める（実測 12 倍） | 実装していない（`xgcd_floor` は素直な形） |
| **§5.4** | **Babai 縮約の実装**。`k` の係数が **30 bit の整数 × 2^s** になるようスケールする。`kf`, `kg` の計算が「縮約の最も高価な部分」。**多倍長×ワードの二次法**と **RNS+NTT** の 2 択があり、閾値は**測って決めよ** | 分岐版の `BABAI_STEP` と `_negacyclic_addmul!`。`docs/debug_log.md` #044 |
| §5.5 | 「支配的なコストは Babai 縮約のままである」 | 実測でも同じ（#040, #044） |

## 到達できていないもの

| 資料 | 状態 |
|:---|:---|
| **FIPS 206 ドラフト**（FN-DSA） | 未入手。実行環境の egress ポリシーが nvlpubs.nist.gov を遮断している（`docs/debug_log.md` #002）。このリポジトリが「FN-DSA」を名乗る根拠は Falcon 仕様書 v1.2 の方であり、FIPS 206 が Falcon に加えた変更（ドメイン分離、`ctx` など）は**反映していない** |
| 提出パッケージの `Supporting_Documentation/additional/` | 未入手。`parameters.py`（パラメータ導出の自動化）と `test-vector-sampler-falcon{512,1024}.txt`（より詳細な SamplerZ ベクタ）が入っているはず。前者があれば σmin の ε が `[derived]` でなくなる |

**#002 の「仕様書 PDF に到達できない」は 2026-08-21 に解消した**（#046）。
`FalconParams.spec_ref` は Table 3.3 を指しており、
`test/test_spec.jl` が表そのものを転記して検査している。

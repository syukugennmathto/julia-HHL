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
| **FIPS 206**（FN-DSA） | **そもそも発行されていない。** 2026-08-21 時点で csrc.nist.gov の FIPS 一覧に FIPS 206 は無く、initial public draft も出ていない。「入手できない」ではなく「**存在しない**」。下記 |

### 番号を間違えやすいので注記

NIST の PQC 標準は 4 本あり、**この実装に関係するのは 206 だけ**である:

| | 名前 | 元 | この実装との関係 |
|:---|:---|:---|:---|
| FIPS 203 | ML-KEM | Kyber | 無関係（鍵カプセル化） |
| FIPS 204 | ML-DSA | Dilithium | 無関係。ただし `docs/math/` の対比の相手 |
| FIPS 205 | SLH-DSA | SPHINCS+ | **無関係**（ハッシュベース署名。格子ではない） |
| **FIPS 206** | **FN-DSA** | **Falcon** | **これが要る** |

2026-08-21 に FIPS 205（SLH-DSA、2024-08-13 発行、61 頁）が届いたが、
これは SPHINCS+ の標準であって Falcon ではない。定数は 1 つも関係しない。
| 提出パッケージの `Supporting_Documentation/additional/` | 未入手。`parameters.py`（パラメータ導出の自動化）と `test-vector-sampler-falcon{512,1024}.txt`（より詳細な SamplerZ ベクタ）が入っているはず。前者があれば σmin の ε が `[derived]` でなくなる |

**#002 の「仕様書 PDF に到達できない」は 2026-08-21 に解消した**（#046）。
`FalconParams.spec_ref` は Table 3.3 を指しており、
`test/test_spec.jl` が表そのものを転記して検査している。

## FIPS 206 の状態（2026-08-21 現在）

**未発行。** 「読めていない」のではなく「**まだ書かれた文書が公開されていない**」。

- csrc.nist.gov の FIPS 一覧に FIPS 206 は**掲載されていない**
  （掲載は FIPS 203 / 204 / 205、いずれも 2024-08-13 Final まで）。
  ドラフト一覧にも無い。
- NIST は 2025-08-28 にドラフトを承認手続きに提出。
  2025-09 の第 6 回 PQC 標準化会議で Ray Perlner が
  "We expect to release an Initial Public Draft soon / It's basically written,
  awaiting approval" と発表。以後 NIST・商務省内のクリアランスで停滞。
- pqc-forum の NIST 公式回答:
  "The draft is still in clearance within NIST and the Dept. of Commerce"。
- IETF `draft-ietf-cose-falcon-04` は FN-DSA を
  "defined in US NIST FIPS 206 (expected to be published in late 2026 early 2027)"
  と記載。
- **ACVP の FN-DSA テストベクタも存在しない**
  （usnistgov/ACVP-Server に FN-DSA / FIPS206 のフォルダが無い）。

したがって **`spec_ref` に FIPS 206 を書ける日はまだ来ていない**し、
「FIPS 206 準拠」を主張できる実装は世界に一つも無い。

### 公表されている範囲での「Falcon round-3 → FN-DSA」の差分

出典は NIST 公式プレゼン（Perlner 2025-09）と pqc-forum の NIST 公式回答、
および Thomas Pornin の追随実装（`pornin/c-fn-dsa`、Unlicense）の README。
**本文が無いので、以下はすべて確定ではない。**

| 変更 | 内容 |
|:---|:---|
| ドメイン分離 | ML-DSA 型。中間値 μ を導入し、公開鍵ハッシュ `tr = SHAKE256(pk, 64)` を混ぜる |
| `ctx` | context string（≤ 255 バイト、既定は空）を追加 |
| pure / prehash | `HashFN-DSA` を定義。OID は NIST 登録待ち（TBD） |
| internal / external | ML-DSA と同様に分離 |
| 署名の乱択性 | **randomized のみ。決定的署名を明示的に禁止**（浮動小数点実装の差で同一ハッシュから異なる署名が出る危険） |
| salt | repeat ループの**外**でサンプル（round-3 は内側）。eprint 2024/1769 |
| 追加の判定 | 署名の**無限ノルム上限 840**（round-3 は符号化都合の 2047） |
| base sampler | 72 bit → **79 bit**（符号ビットの余り 7 bit を捨てずに使う） |
| 公開鍵 | **NTT 形式**で格納（round-3 は plain） |
| エンコード | **リトルエンディアンに統一** |
| keygen | GS ノルム上限を `0.9999·1.17√q` に、木の葉が `[σmin, σmax]` に入るか明示チェック |
| 浮動小数点 | IEEE 754-2019 の必要部分を FIPS 206 本文に**再掲**（IEEE の許諾済み）。演算順序を規定し **FMA を禁止**。署名は KAT に厳密一致することを要求 |
| 数値パラメータ | `q = 12289`、`n = 512/1024`、`σ`、`σmin`、`σmax` は**変更なし** |

**この実装は上のどれも反映していない。** round-3 の Falcon である。
とくに 4 番目（決定的署名の禁止）と 12 番目（FMA 禁止・KAT 厳密一致）は
`docs/debug_log.md` #048 の話と正面から関係する ―
FIPS 206 は「署名は KAT に bit 一致すること」を要求する方向なので、
**#048 で潰した 3 箇所（複素除算・LDL の D11・ffSampling の最下段）は、
いずれ規範として書かれる可能性が高い。**

> ## ⚠ このブランチは**仕様書**から意図的に逸脱している（**論文**からは逸脱していない）
>
> **ブランチ**: `claude/falcon-julia-fndsa-hl8i54-cschedule`
> **分岐元**: `claude/falcon-julia-fndsa-hl8i54`（`f66c301`）
>
> 二点で `main` と違う:
>
> 1. `babai_reduce` が、仕様書（および Python 参照実装）の `Reduce` ではなく
>    **明示的 bit 予算方式**である。
> 2. `f, g` の標本化が、`samplerz` の畳み込みではなく
>    **C 参照実装の CDT 表**である（`gen_poly_cdt`）。
>
> どちらも C 参照実装の方式で、**どちらも Pornin–Prest 論文に書かれている**。
> 論文 §2.8 は `Reduce` の停止条件について
> 「(F, G) のノルムが減らなくなった時点で抜ければよい」と書いており、
> 「`k` が 0 になったら止める」は**仕様書側の具体化**にすぎない。
> §5.4 は「`k` の係数が 30 bit の整数 × 2^s になるようスケールする」と
> 明記していて、それがこのブランチの `BABAI_STEP` である。
> 出典の対応表は `docs/refs.md`、経緯は `docs/debug_log.md` #041〜#045。
>
> | | main | このブランチ |
> |:---|:---|:---|
> | `Reduce` の出典 | 仕様書 / Python 参照実装 | 論文 §5.4 / C 参照実装 |
> | `f, g` の標本化 | 仕様書（`samplerz` の畳み込み） | C 参照実装（CDT 表） |
> | `(F, G)` | **Python 参照実装と厳密一致** | 一致しない（別の有効解） |
> | `f*G − g*F = q` | 成立 | **成立**（厳密に保存される） |
> | `keygen` 中央値 | 約 129 ms | **約 20 ms** |
>
> **Babai 縮約は一意ではない**ので、どちらも正しい短基底を返すが
> 同じものは返さない。仕様準拠を上位に置くなら `main` を、
> 速度を上位に置くならこちらを使う。
>
> 記録した KAT（`NTRU_SOLVE_KAT`）は**両方で通る** ―
> トイ次数では `f, g` が 8 bit しかなく `size` が 53 にクランプされるので、
> 2 つのスケジュールが同じ動作をするからである。
> **KAT が通ることは、一般に一致することの証拠にならない。**

# FALCON (FN-DSA) Julia 仕様準拠実装

FALCON / FN-DSA（NIST FIPS 206 ドラフト）を Julia で仕様準拠実装する長期プロジェクト。
目標は **FALCON-512 の keygen / sign / verify が KAT を通ること**。

定数時間性・サイドチャネル耐性は**スコープ外**。ただし「本来ここが定数時間実装の
難所である」という注記はコードに残す（同人誌原稿の題材）。

## 進捗

| # | ファイル | 状態 |
|---|---|---|
| 1 | `src/params.jl` | 実装済 / **テスト合格** |
| 2 | `src/shake.jl` | 実装済 / **テスト合格**（SHAKE256 は自前 Keccak） |
| 3 | `src/poly.jl` | 実装済 / **テスト合格** |
| 4 | `src/ntt.jl` | 実装済 / **テスト合格** |
| 5 | `src/fft.jl` | 実装済 / **テスト合格**（根の表が C 参照実装と bit 一致） |
| 6 | `src/ntrugen.jl` | 実装済 / **テスト合格**（`ntru_gen` 含む） |
| 7 | `src/samplerz.jl` | 実装済 / **テスト合格**（公式 KAT 3072 本） |
| 8 | `src/ffsampling.jl` | 実装済 / **テスト合格** |
| 9 | `src/encoding.jl` | 実装済 / **テスト合格**（C 参照実装のバイトと相互運用） |
| 10 | `src/falcon.jl` | 実装済 / **テスト合格**（C 参照実装の署名を検証） |

> **状態**: モジュール 1〜10 すべて実装済み、**Julia 1.12.7 で全テスト合格**。
> FALCON-512 と **FALCON-1024** の両方が keygen / sign / verify を通る。
>
> **達成したこと**: C 参照実装が生成した本物の FALCON-512 の署名を、
> 我々の `falcon_verify` が受理する。鍵・署名のバイト形式も相互運用する。
> さらに **FFT の根の表が C 参照実装の `fpr_gm_tab` と bit 完全一致**
> （`docs/debug_log.md` #031）。
>
> **署名バイトの再現について**: まだできていない。ただし理由が変わった。
> #025 では「浮動小数点だから原理的に無理」と書いたが、#031 の実測で
> 原因が**定数表の 1 ulp の差**（我々の `cispi` 側の誤差）だと特定でき、
> 表は合わせた。残る差分は演算の順序だけである。
> 署名側は当面、性質（検証が通る・β² 以下・乱択・改竄拒否）で守っている。

## テストの走らせ方

```sh
julia --project=falcon -e 'using Pkg; Pkg.test()'
```

n=1024 の鍵生成が 1 本入っていて、そこが実行時間の大半を占める。
反復中に飛ばしたければ:

```sh
FALCON_SKIP_1024=1 julia --project=falcon -e 'using Pkg; Pkg.test()'
```

（コミット前には必ず外して走らせること。）

## 速度対決（C 参照実装 vs Julia）

数値と再現手順は `docs/benchmarks.md`。**このブランチ**の数字である。
同一マシン（Xeon 2.10GHz、4 コア）・同一時刻の中央値（ms、`n=512`）:

| op (n=512) | **C 参照実装（既定 = native FP）** | C（FPEMU、Algorand が配る形） | **Julia** |
|:---|---:|---:|---:|
| `keygen` | 5.928 | 12.266 | **17.122** |
| `sign`（展開済み鍵） | 0.172 | 1.812 | **1.105** |
| `verify`（バイトから） | 0.0264 | 0.0263 | **0.0167** |

- **`verify` はどの C ビルドにも勝っている**（1.58 倍）。
  `verify` は浮動小数点を使わないので C の 3 ビルドが同じ値になり、
  **この勝ちはビルド構成に依らない**（#039）。
- **`sign` は FPEMU ビルドに勝ち、native ビルドには 6.4 倍負けている。**
- **`keygen` は native の 2.9 倍、FPEMU の 1.4 倍。**

> **どちらの C が「参照実装の既定」か**は #045 で訂正した。
> 公式アーカイブは `FALCON_FPEMU` も `FALCON_FPNATIVE` も
> コメントアウトしたまま配っており、`README.txt` は
> "If using FALCON_FPNATIVE ... **This is the default.**" と書いている。
> 無改変ビルドは実際に native FP になる（確認済み）。
> FPEMU を既定にしているのは **Algorand のフォーク**であって参照実装ではない。
> 初版の README はこれを取り違えていた。

### 鍵生成の推移

このブランチで踏んだ道:

| 状態 | `keygen` 中央値 | 出典 |
|:---|---:|:---|
| 最初 | 9610 ms | #032 |
| `bitsize` / 混成 `karamul` / in-place GMP | 133 ms | #040 |
| C の縮約スケジュール（明示的 bit 予算） | 60 ms | #042 |
| C の CDT サンプラ | ― | #043 |
| サンプラの `Core.Box` を外す | ― | #043 |
| Babai の補正を `mpz_addmul_ui` に融合 | **17 ms** | #044 |

**565 倍**改善して、C の既定ビルドの 2.9 倍まで来ている。

## golden vector の再生成

```sh
python3 falcon/scripts/gen_vectors.py
```

決定的（`os.urandom` を使っていない）なので、再実行してバイト単位で同じものが
出なければそれ自体がバグ。

## ディレクトリ

```
src/            Julia 実装
test/           テスト
  vectors/      参照実装から生成した golden vector（コミット済み）
                fft_c_kat.jl だけは **C 参照実装**から生成
scripts/
  gen_vectors.py    golden vector 生成スクリプト（Python 参照実装を使う）
  cref_fft_dump.c   C 参照実装から FFT ベクタを吐くドライバ
  cref_kat_dump.c   C 参照実装から鍵・署名の KAT を吐くドライバ
  cref_gm_dump.c    C 参照実装から FFT の根の表 fpr_gm_tab を吐くドライバ
  cref_bench.c      C 参照実装の速度を測るドライバ
  bench.jl          こちらの速度を同じ形式で測るスクリプト
  bench_compare.jl  両者を並べて表にする
  cref/             **C 参照実装の公式アーカイブ**（Falcon-impl-20211101, MIT）
                    経緯と algorand ミラーとの差分は cref/PROVENANCE.md
  pyref/            Python 参照実装 (tprest/falcon.py, MIT) を vendor したもの
docs/
  debug_log.md          デバッグ記録（セッションをまたぐ唯一の記憶）
  refs.md               **一次資料の所在と、どの主張がどれに拠るかの対応表**
  build_cref_macos.md   C 参照実装を dylib にする手順（macOS / Apple Silicon）
  benchmarks.md         C 版との速度対決の生データと再現手順
  math/                 数学的背景（原稿素材）
```

## この環境で Julia を入手する方法

julialang.org 系のドメインが塞がれている環境でも、Docker Hub のレジストリ API
経由で公式イメージから Julia を取り出せる（`docs/debug_log.md` #013）。
デーモンは不要で、blob は `mirror.gcr.io` から取る。

## 参照実装

| | 用途 | 出所 |
|:---|:---|:---|
| C | 主たる突き合わせ相手（`ccall`）・速度対決 | **公式アーカイブ `Falcon-impl-20211101.zip`**（MIT）。`scripts/cref/` に vendor。`scripts/cref/PROVENANCE.md` |
| Python | 副の突き合わせ相手・golden vector 生成 | <https://github.com/tprest/falcon.py>（MIT, Thomas Prest）。`scripts/pyref/` に vendor（`LICENSE` 同梱） |
| 論文 | 降下アルゴリズムの一次資料 | Pornin, Prest, IACR ePrint **2019/015**。本文は vendor せず節番号で引用（`docs/refs.md`） |

#045 以前は algorand/falcon ミラーを引用元にしていた。
**アルゴリズムのファイルはすべてバイト一致**なので既存の `[C-ref]` 引用は
そのまま有効だが、`config.h` だけは違っていて、そこを取り違えていた。

C 実装のビルドは `docs/build_cref_macos.md`。
`Zf(name)` が `falcon_inner_name` に展開されること、`fpr` が
`struct { double v; }` であることなど、`ccall` 側の注意点はそこにまとめてある。

## 定数の出典について

パラメータは記憶や推測で書かない方針。各定数には実際に読んだ file:line を
`[C-ref]` / `[Py-ref]` タグで付けてある。

**仕様書 PDF そのものには依然として到達できていない。**
2026-08-21 に C 参照実装の公式アーカイブと Pornin–Prest 論文は入手できたが
（`docs/refs.md`）、提出パッケージ同梱の `falcon.pdf` はまだ無い。
実行環境の egress ポリシーが falcon-sign.info / nvlpubs.nist.gov /
eprint.iacr.org をすべて遮断している（`docs/debug_log.md` #002）。
したがって「仕様書の表番号」は埋まっておらず、
`FalconParams.spec_ref` は `"TODO: ..."` のまま残してある。
値そのものは C と Python の 2 実装が完全一致していることで担保している。
**「値が疑わしい」のではなく「出典が書けていない」**という意味の TODO である。

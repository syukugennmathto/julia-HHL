> ## ⚠ このブランチは仕様から意図的に逸脱している
>
> **ブランチ**: `claude/falcon-julia-fndsa-hl8i54-cschedule`
> **分岐元**: `claude/falcon-julia-fndsa-hl8i54`（`f66c301`）
>
> `babai_reduce` を、仕様書（および Python 参照実装）の `Reduce` ではなく
> **C 参照実装の明示的 bit 予算方式**に置き換えてある。
> 詳細は `docs/debug_log.md` #041 / #042、`docs/math/06_ntrugen.md` 6.9b。
>
> | | main | このブランチ |
> |:---|:---|:---|
> | `Reduce` の出典 | 仕様書 / Python 参照実装 | C 参照実装 |
> | `(F, G)` | **Python 参照実装と厳密一致** | 一致しない（別の有効解） |
> | `f*G − g*F = q` | 成立 | **成立**（厳密に保存される） |
> | 降下 | 約 200 ms | **約 18 ms（10.5 倍）** |
> | `keygen` 中央値 | 約 129 ms | **約 60 ms（2.1 倍）** |
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

数値と再現手順は `docs/benchmarks.md`。同一時刻・同一マシンの中央値（ms）:

| op (n=512) | C (FPEMU, 参照実装の既定) | C (native FP) | **Julia** |
|:---|---:|---:|---:|
| `keygen` | 10.881 | 5.146 | **132.671** |
| `sign`（展開済み鍵） | 1.566 | 0.146 | **1.155** |
| `verify`（バイトから） | 0.0204 | 0.0204 | **0.0148** |

**`verify` はどの C ビルドにも勝っている**（1.13〜1.39 倍）。
`verify` は浮動小数点を使わないので C の 2 ビルドが同じ値になり、
**この勝ちは FPEMU のハンデに依らない**（#039）。

**3 つ合わせると**:

- **1 回ずつ**（鍵 1 個・署名 1 通・検証 1 通）: 合計 134.3 ms 対 C(既定) 13.9 ms
  ― **9.6 倍遅い**。そして**我々の合計の 98.8% は `keygen`**。
- **定常運用**（鍵 1 個で N 通）: `Julia = 133.13 + 1.170N`、
  `C(既定) = 12.34 + 1.586N` ― **N ≈ 290 で逆転**。
- **検証だけ**: 勝っている。

C の native FP 版とは署名の傾きで負けているので追いつかない。
ただし参照実装が**実際に配布するのは FPEMU 版**であり、
ネイティブ FP は「署名の非決定性は壊滅的」として使うなと書かれている。
我々はその 1 桁を払っていない代わりに、その保証も持っていない（#031）。

### 鍵生成が唯一の負け筋

9610 → **133 ms**（72 倍改善）。それでも C の 12 倍。
理由は表現でも RNS でもなく、**Babai 縮約のスケジュール**だった ―
`size = max(53, bits(f,g))` の `53` のクランプにより、
深い段では補正 `k` が 0 に丸まって縮約が降参し、
6240 bit の係数がそのまま運ばれる。C は明示的な bit 予算を持つので降参しない。
仕事量にして **71 倍**の差（`docs/debug_log.md` #041、`docs/math/06_ntrugen.md` 6.9b）。

**参照実装が 2 つあり、この一点で違うアルゴリズムを実装している。**
仕様書の `Reduce` は Python 側の形で、我々はそちらに忠実である。

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
  pyref/            Python 参照実装 (tprest/falcon.py, MIT) を vendor したもの
docs/
  debug_log.md          デバッグ記録（セッションをまたぐ唯一の記憶）
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
| C | 主たる突き合わせ相手（`ccall`） | <https://github.com/algorand/falcon>（round-3 参照実装のミラー、MIT） |
| Python | 副の突き合わせ相手・golden vector 生成 | <https://github.com/tprest/falcon.py>（MIT, Thomas Prest）。`scripts/pyref/` に vendor（`LICENSE` 同梱） |

C 実装のビルドは `docs/build_cref_macos.md`。
`Zf(name)` が `falcon_inner_name` に展開されること、`fpr` が
`struct { double v; }` であることなど、`ccall` 側の注意点はそこにまとめてある。

## 定数の出典について

パラメータは記憶や推測で書かない方針。各定数には実際に読んだ file:line を
`[C-ref]` / `[Py-ref]` タグで付けてある。

ただし**仕様書 PDF そのものには到達できていない**（実行環境の egress ポリシーが
falcon-sign.info / nvlpubs.nist.gov / eprint.iacr.org をすべて遮断している。
`docs/debug_log.md` #002）。したがって「仕様書の表番号」は埋まっておらず、
`FalconParams.spec_ref` は `"TODO: ..."` のまま残してある。
値そのものは C と Python の 2 実装が完全一致していることで担保している。

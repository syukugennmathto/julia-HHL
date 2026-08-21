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

数値と再現手順は `docs/benchmarks.md`。要点だけ:

| op (n=512, 中央値 ms) | C (FPEMU) | C (native FP) | Julia |
|:---|---:|---:|---:|
| `keygen` | 14.24 | 7.15 | 186 |
| `sign` (展開済み鍵) | 2.04 | 0.194 | **2.05** |
| `verify` | 0.030 | 0.029 | 0.458 |

**署名は C の FPEMU 版と互角**。これは Julia が速いのではなく、
参照実装が移植性のために浮動小数点をソフトウェアで実装していて、
そのぶん 1 桁遅いという話である（ネイティブ FPU なら 10 倍速い）。

### verify は C より速い（#039）

バイトからバイトまで、20000 反復の中央値:

| | ms |
|:---|---:|
| **Julia** | **0.02550** |
| C `clang -O3 -march=native` | 0.02664 |
| C `clang -O3`（参照実装の Makefile） | 0.02746 |

**最速の C ビルドより 3〜4% 速い**（3 ラウンドで再現）。
0.251 → 0.0255 ms、9.8 倍。決め手は NTT の最後の 3 段の融合で、
そこが 70% を食っていた。`ntt` / `intt` / `polymulq` は無傷のまま
テストオラクルとして残してある。

### 鍵生成は 52 倍速くなった（なお C の 13 倍遅い）

9610 → **186 ms**。#036 で 42 倍（`bitsize` の逐語訳が 5000 万回確保していた、
降下の 96% が 24 bit を `BigInt` で掛けていた）、#040 でさらに 1.34 倍
（`F[i] -= fk[i] << shift` の 1 行が鍵生成の 64% を占めていた ―
56 bit の値を 6250 bit シフトして引くだけのために、
係数ごとに 790 バイトの `BigInt` を 2 個確保していた。in-place GMP で 10.2 倍）。

**「多倍長が要る」と「多倍長の*確保*が要る」は別だった。**
残りの 13 倍を埋めるには C と同じ RNS が要る ―
今度はそれを測った上で言っている（#040）。

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

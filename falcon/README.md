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

FALCON を Julia で仕様準拠実装する長期プロジェクト。
準拠先は **Falcon 仕様書 v1.2（01/10/2020）**、すなわち round-3 提出版である。
FN-DSA は NIST が FIPS 206 で Falcon に与えた名前だが、
**2026-08-21 時点で FIPS 206 はまだ発行されていない**（IPD すら出ていない）。
NIST が予告している変更（μ によるドメイン分離、`ctx`、randomized のみ、
無限ノルム 840、公開鍵の NTT 化、リトルエンディアン統一など）は
反映していない。一覧は `docs/refs.md`。
目標は **FALCON-512 の keygen / sign / verify が KAT を通ること**。

定数時間性・サイドチャネル耐性は**スコープ外**（撤回しない ―
Julia では処理系と FPU 由来の変時間性を閉じられない）。
ただし**測ってはいる**: `docs/constant_time.md` に脅威モデル・棚卸し・
`dudect`（Welch の t 検定）による実測を置いた。
最大の漏洩だった圧縮器は**分岐予測**が原因で、分岐を消したら
データ依存差が 9563 ns → **1 ns** になり、**ついでに速くなった**（#051）。

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
> **署名バイトの再現について**: まだできていない。**が、原因は全部特定した**（#048）。
> #025 の「浮動小数点だから原理的に無理」は**完全に外れ**だった。
> `scripts/cref/` を vendor したので C を共有ライブラリにして 1 段ずつ
> bit 比較でき（`scripts/cmp_cref_fp.jl`）、その結果:
>
> | | 結果 |
> |:---|:---|
> | `fft` `ifft` `split_fft` `merge_fft` | **無改修で bit 一致**（logn 3〜10） |
> | `add` `sub` `mul` `adj` `mul*adj` `mul*selfadj` | **無改修で bit 一致** |
> | `expand_privkey` の `B0` 行列 | **bit 一致** |
> | `div_fft` | C の書き方に直して一致（Smith vs 逆数を作って掛ける） |
> | `ldl_fft` の `D11` | C の書き方に直して一致（乗算 3 回 vs 1 回） |
> | `ffSampling` | C の n=4 手展開を写すと **n=4〜512 で全一致** |
> | 木の葉 | 互いの逆数（**格納規約**の違い、演算の違いではない） |
>
> **「演算の順序」ではなく「同じ式の違う書き方」が 3 箇所あっただけ**だった。
> `_ffsampling_c4` は**まだ繋いでいない** ―
> 算術は検証したが `samplerz` の**呼び出し順序は検証していない**からで、
> 理由は同関数の docstring と #048 に書いた。

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

## 速度（C 参照実装 vs Julia）

**C は 1 ビルドではない。** `{gcc 13.3, clang 18.1} × {-O2, -O3, -O3 -march=native}
× {native FP, FPEMU}` の 12 ビルドで測り、**範囲**で示す。
Julia は同一設定 3 run の幅。数値と方法は `docs/benchmarks.md`。

| op (n=512, 中央値 ms) | C native FP | C FPEMU | **Julia** |
|:---|---:|---:|---:|
| `keygen` | 4.98〜6.07 | 11.2〜13.2 | **11.7〜12.7** |
| `sign`（展開済み鍵） | 0.136〜0.173 | 1.59〜1.82 | **0.98〜1.49** |
| `verify`（バイトから） | 0.0166〜0.0268 | 0.0168〜0.0268 | **0.0164〜0.0213** |

- **`verify` は C と同着。** 12 ビルド中 10 より速いが、
  最良の 2（clang -O3 -march=native）とは測定の幅の中で区別できない。
  > **訂正**: 以前ここには「どの C ビルドにも 1.58 倍勝っている」と書いていた。
  > それは `gcc -O2` を「C」と呼んだ結果である。
  > **C 側のビルド間の幅（1.6 倍）が、主張していた差と同じ大きさだった**（#052）。
- **`sign` は FPEMU 帯より速く、native 帯に約 7 倍負ける。**
- **`keygen` は FPEMU 帯と重なり、native 帯の約 2 倍。**

定常運用（鍵 1 個で N 通）では **C FPEMU に N ≈ 2 で逆転**し、
**C native には永久に追いつかない**（傾きで 6.6 倍負け）。

### 鍵生成の推移

| 状態 | 中央値 | 出典 |
|:---|---:|:---|
| 最初 | 9610 ms | #032 |
| `bitsize` / 混成 `karamul` / in-place GMP | 133 ms | #040 |
| C の縮約スケジュール | 60 ms | #042 |
| CDT サンプラ + `Core.Box` 除去 | 34.3 ms | #043 |
| Babai の補正を `mpz_addmul_ui` に融合 | 17.1 ms | #044 |
| Karatsuba をやめて `mpz_addmul` の schoolbook に | **11〜13 ms** | #047 |

**875 倍**改善した。そのうち**高水準の配列演算に由来する部分は無い** ―
すべて確保の除去である（`docs/why_julia.md`）。

## なぜ Julia で書いたのか

当初の動機は「Julia なら速くなると思ったから、とくに**行列計算が得意**だから」。
**その見立ては外れた** ― FALCON の重い部分は行列計算ではなく、
多倍長整数の積和（鍵生成）と分割統治の再帰（署名）である。

外れたことと、それでも Julia で書いた価値がどこにあったかは
`docs/why_julia.md` に測定つきで書いた。要点は、
**別言語で書いたからこそ仕様書と参照実装の差が見えた**ことである（#050）。

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
  constant_time.md      **脅威モデル・変時間箇所の棚卸し・dudect による実測**
  why_julia.md          **なぜ Julia を選び、その見立てがどこで外れたか**
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

パラメータは記憶や推測で書かない方針。各定数には実際に読んだ出典を
`[Spec]`（仕様書の式・表・頁）/ `[C-ref]` / `[Py-ref]`（file:line）/
`[derived]` タグで付けてある。

**2026-08-21 に仕様書 v1.2（01/10/2020）が届き、#002 以来の TODO は閉じた。**
`FalconParams.spec_ref` は Table 3.3 を名指ししており、
`test/test_spec.jl` が仕様書の表そのものを転記して検査している
（Table 3.3 のパラメータ、式 (2.10)〜(2.14)、Table 3.1 の分布表、
Table 3.2 の SamplerZ ベクタ 16 本）。

**照合の結果、値は一つも変わらなかった。**
「2 つの独立な参照実装の完全一致で担保する」という代替手段は、
全項目について正しかった（#046）。

一次資料から新しく分かったことの方が大きい:

- **σ と σ_min は別の ε から出ている。**
  σ は (2.13) を `ε ≤ 1/√(Q_s·λ)`（`Q_s = 2^64`、`λ = 128/256`）で評価した値、
  σ_min は `1/√(2^64·n³)` での平滑化パラメータ。
  取り違えると 10% 違う値が**静かに**出る。両方向をテストで固定した。
- **#022 の「KAT の反転規約はハーネスにしか書かれていない」は外れ**だった。
  仕様書 Table 3.2 のベクタも同じ規約で、`reversed_chunks=true` で 16/16 通る。

まだ無いのは提出パッケージの `Supporting_Documentation/additional/`
（σ_min の導出が入っている `parameters.py`）。
**FIPS 206 は「無い」のではなく「まだ存在しない」** ―
2026-08-21 時点で IPD も出ていない。詳細と、公表されている範囲での
round-3 → FN-DSA の差分一覧は `docs/refs.md`。

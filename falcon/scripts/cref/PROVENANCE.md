# FALCON 参照実装（C）― vendor 済み

## 出所

**Falcon reference implementation, `Falcon-impl-20211101.zip`**
（<https://falcon-sign.info/> の "Falcon reference implementation: Source code
archive [zip]"）。ユーザが本セッションに直接アップロードしたもの。
アーカイブ名は 2021-11-01 だが、`README.txt` の `Version:` は `2020-09-30`
で、ファイルの mtime も 2020-10-07 である。2021-11-01 の差し替えは
`shake256_init_prng_from_seed()` / `shake256_init_prng_from_system()` の
**外部 API のバグ修正**（David Lazar, Chris Peikert が発見）によるもので、
NIST 提出物・KAT・ベンチマークは影響を受けない、と公式サイトに注記がある。

ライセンスは **MIT**（各 `.c` の冒頭、`Copyright (c) 2017-2019 Falcon Project`、
`@author Thomas Pornin`）。

## それまで使っていたミラーとの関係

セッション #001〜#043 の `[C-ref]` 引用は
<https://github.com/algorand/falcon>（round-3 参照実装のミラー）を読んで書いた。
公式アーカイブが手に入ったので**全ファイルを差分にかけた**:

| ファイル | 差分 |
|:---|:---|
| `keygen.c` `fft.c` `fpr.c` `sign.c` `codec.c` `common.c` `rng.c` `shake.c` `fpr.h` | **0 行**（バイト一致） |
| `vrfy.c` | 13 行 ― Algorand が `mq_NTT` 等を `deterministic.h` 用に export しただけ |
| `inner.h` | 31 行 ― 同上の宣言 |
| `config.h` | 143 行 ― **後述。ここだけは意味のある違いがある** |

**アルゴリズムを含むファイルはすべてバイト一致**なので、
これまでに書いた `[C-ref] keygen.c:1234` の類の引用は
**公式アーカイブに対してもそのまま有効**である。

## config.h の違い ― これは記録すべき誤解だった

`config.h` で Algorand は `FALCON_FPEMU 1` を**有効にしている**。
理由は彼らの `config.h` に書いてある: 決定的署名を得るためであり、
非決定的な署名は「CATASTROPHIC SECURITY FAILURE」に繋がりうる、と。

**この文章は Algorand のフォークのものであって、Falcon 参照実装のものではない。**

公式アーカイブの `config.h` は `FALCON_FPEMU` も `FALCON_FPNATIVE` も
**両方コメントアウトされている**。そして `README.txt` にはこう書いてある:

> If using FALCON_FPNATIVE, then the C 'double' type is used for all
> floating-point operations. **This is the default.**

実際に公式ツリーを無改変でビルドすると
`# cref_bench fpemu=0 fpnative=1` と出る（確認済み）。

したがって:

- **参照実装が既定でビルドするのはネイティブ FP 版**である。
- `docs/benchmarks.md` と `README.md` の初版が
  「参照実装が実際に配布するのは FPEMU 版であり、ネイティブ FP は
  使うなと書かれている」と書いていたのは**誤り**で、
  それは Algorand のフォークの方針を参照実装の方針と取り違えていた。
  `docs/debug_log.md` #045 に記録した。
- 速度対決の「フェアな相手」は **native FP ビルド**の方である。

## ビルド

```sh
cc -O2 -o cref_bench ../cref_bench.c \
   codec.c common.c falcon.c fft.c fpr.c keygen.c rng.c shake.c sign.c vrfy.c -lm
```

`config.h` を触らなければネイティブ FP 版になる。FPEMU 版が欲しければ
`config.h` の `#define FALCON_FPEMU 1` のコメントを外す
（`-DFALCON_FPEMU=1` はコマンドラインからでは効かない ―
`config.h` 側の `#define` が後から上書きするため。
Algorand ツリーで `-DFALCON_FPNATIVE=1` を渡すと
`warning: "FALCON_FPNATIVE" redefined` が出て**無視される**）。

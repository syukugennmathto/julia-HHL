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
| 5 | `src/fft.jl` | 実装済 / **テスト合格**（C 参照実装とも突き合わせ済み） |
| 6 | `src/ntrugen.jl` | 実装済 / **テスト合格**（`ntru_gen` は module 7 待ち） |
| 7 | `src/samplerz.jl` | 未着手 |
| 8 | `src/ffsampling.jl` | 未着手 |
| 9 | `src/encoding.jl` | 未着手 |
| 10 | `src/falcon.jl` | 未着手 |

> **状態**: モジュール 1〜5 について **9672 件のテストが全て合格**している。
> ただし実行に使ったのは **Julia 1.11.9**（Docker イメージから取り出したもの。
> 経緯は `docs/debug_log.md` #013）であって、**ユーザ環境の 1.12 ではない**。
> 1.12 でも一度確認すること。
>
> 最初の実行では 16 failed / 21 errored だった。内訳と原因は
> `docs/debug_log.md` #014〜#019 に全て記録してある。うち 1 件
> （ChaCha20 のカウンタ溢れ）は**参照実装のほうが間違っていた**。

## テストの走らせ方

```sh
julia --project=falcon -e 'using Pkg; Pkg.test()'
```

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
  pyref/            Python 参照実装 (tprest/falcon.py, MIT) を vendor したもの
docs/
  debug_log.md          デバッグ記録（セッションをまたぐ唯一の記憶）
  build_cref_macos.md   C 参照実装を dylib にする手順（macOS / Apple Silicon）
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

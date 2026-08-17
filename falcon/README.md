# FALCON (FN-DSA) Julia 仕様準拠実装

FALCON / FN-DSA（NIST FIPS 206 ドラフト）を Julia で仕様準拠実装する長期プロジェクト。
目標は **FALCON-512 の keygen / sign / verify が KAT を通ること**。

定数時間性・サイドチャネル耐性は**スコープ外**。ただし「本来ここが定数時間実装の
難所である」という注記はコードに残す（同人誌原稿の題材）。

## 進捗

| # | ファイル | 状態 |
|---|---|---|
| 1 | `src/params.jl` | 実装済 / **Julia 未実行** |
| 2 | `src/shake.jl` | 実装済 / **Julia 未実行** |
| 3 | `src/poly.jl` | 実装済 / **Julia 未実行** |
| 4 | `src/ntt.jl` | 未着手 |
| 5 | `src/fft.jl` | 未着手 |
| 6 | `src/ntrugen.jl` | 未着手 |
| 7 | `src/samplerz.jl` | 未着手 |
| 8 | `src/ffsampling.jl` | 未着手 |
| 9 | `src/encoding.jl` | 未着手 |
| 10 | `src/falcon.jl` | 未着手 |

> **重要**: モジュール 1〜3 を書いたセッションの実行環境には Julia が入って
> おらず（`docs/debug_log.md` #001 参照）、**テストは一度も実行されていない**。
> 期待値（`test/vectors/*.jl`）は参照実装を実際に走らせて生成済みなので信頼できるが、
> Julia コード側は構文チェックすら通っていない。最初にやることは下記のテスト実行。

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
scripts/
  gen_vectors.py  golden vector 生成スクリプト
  pyref/          Python 参照実装 (tprest/falcon.py, MIT) を vendor したもの
docs/
  debug_log.md          デバッグ記録（セッションをまたぐ唯一の記憶）
  build_cref_macos.md   C 参照実装を dylib にする手順（macOS / Apple Silicon）
  math/                 数学的背景（原稿素材）
```

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

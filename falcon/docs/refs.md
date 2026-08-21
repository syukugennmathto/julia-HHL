# 一次資料の所在

このプロジェクトの方針は「パラメータや定数を記憶や推測で書かず、
必ず出典を書く」である。到達できているもの・いないものをここに一覧する。

## 到達できているもの

| # | 資料 | 所在 | 用途 |
|:--|:---|:---|:---|
| A | **Falcon 参照実装（C）** `Falcon-impl-20211101.zip` | `scripts/cref/`（vendor 済、MIT）。経緯は `scripts/cref/PROVENANCE.md` | `[C-ref]` タグの引用元。速度対決の相手 |
| B | **Python 参照実装** tprest/falcon.py | `scripts/pyref/`（vendor 済、MIT） | `[Py-ref]` タグの引用元。golden vector 生成 |
| C | **Pornin, Prest: More Efficient Algorithms for the NTRU Key Generation using the Field Norm** (IACR ePrint **2019/015**, 31 頁) | 本文はリポジトリに含めない（下記） | 降下アルゴリズムの一次資料 |

### C（論文）をリポジトリに入れていない理由

MIT ではないので、全文を vendor しない。
**節番号で引用し、必要な箇所だけ短く引く**という形にしてある
（この文書と `docs/debug_log.md`、`docs/math/06_ntrugen.md` がそれ）。
入手先: <https://eprint.iacr.org/2019/015>。

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
| **Falcon 仕様書 PDF**（`falcon.pdf`、提出パッケージ同梱） | **未入手。** 実行環境の egress ポリシーが falcon-sign.info / nvlpubs.nist.gov / eprint.iacr.org を遮断している（`docs/debug_log.md` #002）。ユーザから届いたのは**参照実装のアーカイブ**であって提出パッケージではなかった |
| **FIPS 206 ドラフト** | 未入手（同上） |

したがって `FalconParams.spec_ref` の「仕様書の表番号」は**まだ埋まっていない**。
値そのものは C と Python の 2 実装が完全一致していることで担保している。
**「値が疑わしい」のではなく「出典が書けていない」**という意味の TODO である。

必要なのは提出パッケージ（`falcon-round3.zip` など）に入っている
`falcon.pdf` の 1 ファイルだけ。

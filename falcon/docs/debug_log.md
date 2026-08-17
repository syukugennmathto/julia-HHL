# FALCON Julia 実装 デバッグ記録

セッションをまたぐ唯一の記憶。**ここに書いていないことは存在しなかったことになる。**

エントリ番号は通し番号。フォーマット:

```markdown
## [YYYY-MM-DD] #NNN module: 一行症状
- **モジュール**:
- **症状**:
- **再現条件**:
- **切り分け**:
- **外した仮説**:      ← 必ず埋める。何も外していないなら「なし」
- **原因**:            ← 説明できないなら「不明（要調査）」と書いて先へ進まない
- **修正**:
- **回帰テスト**:
- **学び**:
```

---

## [2026-08-17] #001 環境: このセッションでは Julia が実行できない

- **モジュール**: （環境）
- **症状**: モジュール 1〜3 を書いたが、`julia --project=falcon -e 'using Pkg; Pkg.test()'`
  を 1 度も実行できていない。つまり**現時点のコードは構文チェックすら通っていない**。
- **再現条件**: リモート実行環境（Ubuntu 24.04 コンテナ）で常に。
- **切り分け**:
  - [x] `which julia` → 無し
  - [x] julialang-s3.julialang.org からの取得 → egress proxy が 403（CONNECT 拒否）
  - [x] GitHub Releases (`JuliaLang/julia`) → Linux バイナリの release asset が無く 404
  - [x] `apt-get install julia` → noble に julia パッケージが存在しない（candidate: none）
  - [x] `cache.julialang.org` / `install.julialang.org` → いずれも proxy が 403
  - [→] 到達可能なのは `raw.githubusercontent.com`, `github.com`,
        `archive.ubuntu.com`, PyPI のみ
- **外した仮説**: 「Ubuntu の universe に julia がいるはず」と思って
  `apt-cache policy julia` まで行ったが、Julia は Debian/Ubuntu アーカイブから
  削除済みで、universe を有効にしても入らない。ここで 10 分溶かした。
- **原因**: 実行環境のネットワークポリシーが julialang.org 系ドメインを許可して
  いない。Julia の公式 Linux バイナリは GitHub Releases では配布されていないため、
  許可されているホストだけでは入手経路が無い。
- **修正**: コードの検証責任を分割した。
  1. **数学の正しさ**はこのコンテナ内で検証する。Python 参照実装を
     `falcon/scripts/pyref/` に vendor し、`falcon/scripts/gen_vectors.py` で
     golden vector を生成して `falcon/test/vectors/*.jl` にコミットした。
     生成時に参照実装を実際に走らせているので、期待値そのものは検証済み。
  2. **Julia コードの正しさ**は macOS 側で `Pkg.test()` を走らせて確認する。
     これは未実施。
- **回帰テスト**: なし（環境の問題であってコードの問題ではない）。ただし
  `test/vectors/*.jl` が生成物としてコミットされているので、Julia が無い環境でも
  テストは再現可能になった。
- **学び**: 「参照実装と突き合わせる」方法論は、**参照実装を動かせる場所**と
  **移植先を動かせる場所**が同じであることを暗黙に仮定していた。分かれている場合は、
  境界に置くのは「実行」ではなく「コミットされたベクタ」にするしかない。
  次のセッションで最初にやることは、macOS 上で `Pkg.test()` を走らせて
  **モジュール 1〜3 の赤を全部拾うこと**。ここで出る失敗は #002 以降になる。

---

## [2026-08-17] #002 params: 仕様書 PDF に到達できず、表番号を引けなかった

- **モジュール**: params.jl
- **症状**: 「パラメータの具体値は必ず仕様書の表から引き、出典（表番号）を
  コメントに書く」という要件を、**表番号の部分だけ満たせていない**。
- **再現条件**: 常に。
- **切り分け**:
  - [x] `falcon-sign.info` → egress proxy が 403
  - [x] `nvlpubs.nist.gov`（FIPS 206 ipd PDF）→ 403
  - [x] `csrc.nist.gov` → 403
  - [x] `eprint.iacr.org` → 403
  - [x] `www.di.ens.fr/~prest/Publications/falcon.pdf` → 403
  - [x] `datatracker.ietf.org`（draft-ietf-cose-falcon）→ 403
  - [→] `raw.githubusercontent.com` は到達可能 → **参照実装のソースは読める**
- **外した仮説**: 「WebFetch ツールなら proxy とは別の経路を通るのでは」と考えて
  試したが、同じ egress proxy を通っており同じく 403。Web 検索は通るのに
  Web 取得は通らないので、一瞬「取得できているのに解析に失敗している」のかと
  勘違いした。ツールのエラーメッセージ（`EGRESS_BLOCKED`）をちゃんと読めば
  すぐ分かった話。
- **原因**: 実行環境のネットワークポリシー。仕様書 PDF のホストが全滅している。
- **修正**: **表番号を推測で書かない**方針を採った。具体的には:
  - 値そのものは、**相互に独立な 2 つの参照実装**から引いた。
    C 参照実装（`algorand/falcon`、round-3 のミラー）と
    Python 参照実装（`tprest/falcon.py`）で、`sig_bound` と `sig_bytelen` が
    完全一致することを確認済み（`34034726` / `666`）。
  - 各定数に `[C-ref] common.c:250` のような **実際に読んだ file:line** を付けた。
  - `FalconParams.spec_ref` フィールドは `"TODO: ..."` のまま残し、
    `test_params.jl` がそれを検査している（`@test occursin("TODO", p.spec_ref)`）。
    PDF が読める環境で埋めたら、このテストが「埋め忘れ」ではなく
    「埋めたのにテストを直し忘れ」で赤くなる。それでよい。
- **回帰テスト**: `test/test_params.jl` の以下:
  - byte 長を C のマクロ（Julia に転記した `c_pubkey_size` 等）と logn=2..10 で照合
  - `sigma_min` が平滑化パラメータ `eta_eps(Z)` と一致（後述 #003）
  - `sigma == 1.17 * sqrt(q) * sigma_min`
  - `sig_bound / (2n sigma^2)` が 1〜2 の間（桁を間違えたら落ちる）
- **学び**: 「出典を書け」という指示の本質は表番号そのものではなく
  **「記憶から書くな」**。到達できる一次資料が変われば出典の形式も変わる。
  ただし**形式を勝手に緩めて黙っていてはいけない**ので、`spec_ref` を
  構造体のフィールドとして残し、テストで見えるようにした。

---

## [2026-08-17] #003 params: sigma_min の閉じた式を逆算で特定した

- **モジュール**: params.jl
- **症状**: バグではない。定数を「記憶から書かない」ために、
  `sigma_min = 1.2778336969128337` が**どの式から来ているのか**を
  仕様書なしで確定させる必要があった。
- **再現条件**: —
- **切り分け**:
  - [x] 平滑化パラメータの標準形 `eta_eps(Z) = (1/pi) sqrt(ln(2(1+1/eps))/2)` を仮定
  - [x] 参照値を代入して `eps` を逆算 → `log2(1/eps) = 45.50000000002967` (n=512),
        `47.00000000002099` (n=1024)
  - [→] 45.5 と 47.0 という**ほぼ厳密な値**が出た。差 1.5 = 1.5*(logn の差 1)、
        切片 32 → `1/eps = 2^32 * n^(3/2)`、すなわち `eps = 1/sqrt(2^64 * n^3)`
- **外した仮説**:
  - 最初 `sigma = 1.17 * sqrt(q/(2n))` だと思い込んでいた（`ntrugen.py` の
    コメントに出てくる式なので）。n=512 で 4.05 にしかならず、参照値 165.7 と
    2 桁違う。**これは f, g のサンプリング幅 `sigma_fg` であって
    署名の幅 `sigma` ではない**。FALCON には σ が 4 種類（`sigma`, `sigma_min`,
    `sigma_max`, `sigma_fg`）あり、しかも `sigma_fg` は「4096 個サンプルして
    `4096/n` 個ずつ畳む」という実装都合の基底幅 `1.43300980528773` を持つ。
    ここを混ぜたのが 1 つ目の外し。
  - 次に `eta` の中の `2` を `2n`（`Z^{2n}` の平滑化）だと思って試したが、
    `log2(1/eps)` が 35.5 / 35.0 と汚い値になった。`Z^1` の式が正しい。
- **原因**: （バグではない）
- **修正**: `falcon_eps(n) = 1/sqrt(2^64 * n^3)` と `smoothing_eta(eps)` を
  params.jl に書き、`sigma = 1.17 * sqrt(q) * sigma_min` も確認した
  （1.17 は `ntrugen.py:232` の Gram-Schmidt 上限 `(1.17^2)*q` と同じ定数）。
  閉じた式は Float64 評価で参照値と相対 3e-13 まで一致する。**完全一致では
  ないので、権威ある値は引き続き参照実装のリテラル**とし、閉じた式は
  「桁を打ち間違えたら落ちる」ためのテストに使う。
- **回帰テスト**: `test/test_params.jl` の
  `"sigma_min is the smoothing parameter of Z"` と
  `"sigma = 1.17 * sqrt(q) * sigma_min"` と `"sigma_fg and its base"`。
- **学び**: σ が 4 つある。**どの Gaussian の話をしているかを毎回言え。**
  格子上の Gaussian（`sigma`）、Z 上の各座標の Gaussian（`sigma_min`〜`sigma_max`）、
  鍵生成の f,g の Gaussian（`sigma_fg`）、そしてその実装上の基底幅
  （`SIGMA_FG_BASE`）は全部別物で、単位も意味も違う。
  `1.17` だけが 3 つをまたいで現れる（Gram-Schmidt ノルムの上限）ので、
  この定数を見たら「どの σ の話か」を確認する合図にする。

---

## [2026-08-17] #004 shake: SHA.jl の SHAKE256 は「squeeze の続き」ができない

- **モジュール**: shake.jl
- **症状**: バグではなく設計上の制約。FALCON は SHAKE256 を XOF として使い、
  hash-to-point では**必要な分だけ吸い出しては棄却する**ので、
  「途中まで squeeze して、足りなければ続きを squeeze する」API が要る。
  `SHA.jl` にはそれが無い（`shake256(data, d)` の一発形しかなく、
  `digest!` は呼ぶたびに padding を再適用するので同じ context を 2 度
  squeeze できない）。
- **再現条件**: —
- **切り分け**:
  - [x] `SHA.jl` の `src/shake.jl` を読んだ。`digest!(ctx, d, p)` は
        `d > blocklen` のときは内部で継続 squeeze するが、
        呼び出しをまたいだ継続はできない構造（`used` フラグと再 padding）。
  - [→] XOF の定義から、`SHAKE256(m, d1)` は `SHAKE256(m, d2)` (d1<d2) の
        **prefix** になる。ならばメッセージを保持して長さを倍々にしながら
        取り直せばよい。
- **外した仮説**: なし（最初からこの方針）。
- **原因**: （バグではない）
- **修正**: `SHAKE256XOF` を「吸収済みメッセージ + 生成済みバッファ + 位置」で
  持ち、足りなくなったら倍の長さで `SHA.shake256` を取り直す実装にした。
  償却 O(N) だが、取り直しのたびにメッセージを再吸収するので厳密には
  O(N + |m| log N)。FALCON のサイズでは無視できる。
- **回帰テスト**: `test/test_shake.jl` の `"SHAKE256 prefix property"`。
  **これは FIPS 202 の性質ではなく `SHA.jl` の実装への仮定**なので、
  hashlib 由来のベクタで明示的に検証する。ここが落ちたら
  `squeeze!` を細工するのではなく Keccak-f[1600] を自前で書くこと、
  とテストのコメントに書いた。
- **学び**: 標準ライブラリで代用するとき、**代用が壊れる条件をテストに書いておく**。
  「たぶん prefix になっているはず」で通してしまうと、後で hash-to-point が
  ときどきずれるという最悪の形で返ってくる。

---

## [2026-08-17] #005 shake: 参照 PRNG の randombytes は「余りを捨てる」

- **モジュール**: shake.jl
- **症状**: バグではなく、移植前に潰しておいた罠。
  `scripts/pyref/rng.py` の `randombytes(k)` は、バッファ残量が `k` に
  満たないとき**残りを連結せず捨てて**新しい 512 バイトを生成する。
- **再現条件**: —
- **切り分け**:
  - [x] `rng.py:111-122` を精読
  - [x] `out = "".join(out[i:i+2] for i in range(2*k-2,-1,-2))` のあと
        `bytes.fromhex(out)[::-1]` している。**この 2 つの反転は打ち消し合って
        恒等写像**。実質「バッファ先頭 k バイトを返して k 進める」だけ。
  - [→] 意味があるのは反転ではなく `if 2*k > len(self.hexbytes): 再生成` の分岐
- **外した仮説**: 反転処理に意味があると思って 15 分読んだ。
  C 参照実装がワード単位で消費するのに合わせた名残で、
  バイト列としては何もしていない。**「変な処理があったら、まず恒等でないか
  確かめる」**べきだった。
- **原因**: （バグではない）
- **修正**: Julia 側は素直に「バッファ先頭 k バイト、足りなければ捨てて再生成」で
  実装し、反転は再現しない。その主張は KAT で担保する。
- **回帰テスト**: `test/vectors/chacha20_kat.jl` に
  `[511, 8]` と `[255,255,255]` という**境界をまたぐ要求パターン**を入れた。
  加えて `test_shake.jl` の `"ChaCha20 buffer boundary behaviour"` が
  「捨てている」ことを直接主張している。
- **学び**: 参照実装の「意味の無いコード」と「意味のあるコード」を分けてから移植する。
  ただし**分けた根拠はテストにする**。恒等だと判断した処理が実は恒等でなかった
  場合、それは KAT でしか気づけない。

---

## 次のセッションでやること

1. macOS 上で `julia --project=falcon -e 'using Pkg; Pkg.test()'` を走らせる。
   **モジュール 1〜3 は一度も実行されていない**ので、ここで構文エラー・
   型エラー・テストの赤が出るのが正常。出た分は #006 以降に記録する。
2. C 参照実装を `docs/build_cref_macos.md` の手順で `libfalcon.dylib` にする。
3. モジュール 4（`ntt.jl`）へ進む。

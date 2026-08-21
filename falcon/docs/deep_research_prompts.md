# Deep Research 用プロンプト（別チャットに貼る用）

このリポジトリで**不足している一次資料 2 件**を取りに行くためのプロンプト。
文脈を知らないチャットに貼っても成立するように自己完結させてある。

---

## プロンプト A — FIPS 206 / FN-DSA（優先度: 高）

```
NIST の FIPS 206（FN-DSA: Fast-Fourier Transform over NTRU-Lattice-Based
Digital Signature Algorithm）について、一次資料に基づいて調べてください。
FN-DSA は Falcon を NIST が標準化したものです。

## 私が既に持っているもの（これらは調べ直さなくてよい）

- Falcon 仕様書 v1.2（01/10/2020）— round-3 提出版の specification PDF、67 頁
- Falcon 参照実装 C（Falcon-impl-20211101、MIT）
- Pornin & Prest, "More Efficient Algorithms for the NTRU Key Generation
  using the Field Norm"（IACR ePrint 2019/015）

## 知りたいこと

1. **FIPS 206 の現在の状態**（initial public draft は出ているか、出ているなら
   発行日・版・文書番号は何か、final はいつ予定/発行されたか）。
   NIST の公式ページ（csrc.nist.gov, nvlpubs.nist.gov）の URL を明示すること。

2. **FIPS 206 が Falcon round-3 提出版から変えた点の一覧**。特に:
   - **ドメイン分離** — 署名対象の作り方。round-3 の Falcon は
     `HashToPoint(salt || message)` だが、ML-DSA（FIPS 204）は
     `M' = 0x00 || len(ctx) || ctx || M` のような前置を入れる。
     FN-DSA は同じことをするのか、するならバイト列の正確な定義は何か。
   - **`ctx` 文字列**（context string）の有無・長さ制限・既定値
   - **プレハッシュ版**（HashFN-DSA / pre-hash variant）の有無と、その OID
   - 「内部関数 / 外部関数」（internal vs external）の分離があるか
   - パラメータの改名（Falcon-512 → FN-DSA-512 等）と、
     **数値パラメータそのものに変更があるか**（σ, σ_min, σ_max, ⌊β²⌋,
     q = 12289, 鍵長 897/1793, 署名長 666/1280）。
     **変わっていないなら「変わっていない」と明言してほしい。**
   - 鍵・署名のエンコード形式の変更（ヘッダバイト、salt 長 40、
     Golomb-Rice 圧縮のパラメータ）
   - 決定的署名 / 非決定的署名の扱い、および浮動小数点の再現性について
     FIPS 206 が何を要求しているか（round-3 の参照実装は既定でネイティブ
     double を使い、FPEMU はオプションだった）

3. **FIPS 206 の KAT / ACVP テストベクタ**がどこで公開されているか。
   NIST ACVP (github.com/usnistgov/ACVP-Server) に FN-DSA の
   test vector があるなら、そのパスと形式。

4. **既存の FN-DSA 実装**で FIPS 206 ドラフトに追随しているものがあるか
   （PQClean, liboqs, BoringSSL, pq-crystals, Thomas Pornin の新しい実装など）。
   あるならリポジトリと、Falcon round-3 との差分が読めるコミット/ファイル。

## 出力の形式

- **一次資料（NIST の PDF / 公式リポジトリ）を最優先**。ブログや二次解説は
  一次資料への導線としてのみ使い、事実の根拠にはしないこと。
- 各主張に URL と、可能なら**節番号・アルゴリズム番号・頁**を付けること。
- 「見つからなかった」ことも明記してほしい。曖昧にぼかさないこと。
- 最後に、**「Falcon round-3 に忠実な実装を FIPS 206 準拠にするために
  何を変える必要があるか」のチェックリスト**にまとめてください。
```

---

## プロンプト B — Falcon 提出パッケージの補助資料（優先度: 低）

```
Falcon（post-quantum signature, NIST PQC round 3）の **提出パッケージ**に
同梱されている補助資料の中身について調べてください。

## 探しているファイル

Falcon 公式サイト <https://falcon-sign.info/> が配布している
"Falcon submission package [zip]"（specification, source code, scripts and
test vectors を含む）の中の

    Supporting_Documentation/additional/

というディレクトリ。特に次の 2 つ:

1. **`parameters.py`** — Falcon 仕様書 v1.2 の §2.6 が
   「The resulting parameter selection process is automatized in
   Supporting_Documentation/additional/parameters.py, which also gives the
   core-SVP hardness of key recovery and forgery」と名指ししているスクリプト。

2. **`test-vector-sampler-falcon512.txt` / `test-vector-sampler-falcon1024.txt`**
   — 仕様書 §3.9.3（p.44）が「this submission package contains more extensive
   and detailed test vectors」として名指ししている SamplerZ のテストベクタ。

## 知りたいこと

1. これらのファイルが**オンラインで読める場所**（GitHub のミラー、
   web.archive.org のスナップショット、NIST の round-3 submission
   アーカイブ <https://csrc.nist.gov/projects/post-quantum-cryptography/
   post-quantum-cryptography-standardization/round-3-submissions> など）。
   直リンクを示してほしい。

2. **`parameters.py` が σ_min をどう計算しているか。**
   これが本題です。Falcon 仕様書 Table 3.3 は
   σ_min = 1.277833697（n=512）/ 1.298280334（n=1024）という**値だけ**を
   載せていて、導出式も、使った ε も印刷していない。

   私は逆算して、Z の平滑化パラメータ

       eta_eps(Z) = (1/pi) * sqrt( (1/2) * ln(2 * (1 + 1/eps)) )

   を **eps = 1 / sqrt(2^64 * n^3)** で評価すると相対誤差 3e-13 で
   両方の値が再現することを確認しました。
   （なお仕様書の式 (2.13) が**署名の**標準偏差 σ に使う ε は
   `eps <= 1/sqrt(Qs * lambda)`, Qs = 2^64, lambda = 128/256 で、
   これは σ_min の ε とは**別物**です。式 (2.13) の ε で
   eta を評価すると σ_min から 11% ずれます。）

   知りたいのは、**`parameters.py` が実際にどの式・どの ε を使っているか**。
   私の逆算が当たっているか、それとも別の導出なのか。
   ソースの該当行を引用してほしい。

3. `test-vector-sampler-falcon{512,1024}.txt` の**形式**
   （1 行に何が入っているか: μ, σ', randombytes の hex, 出力 z の順か、
   何本あるか、バイト列を食わせる順序の規約はどうか）。

   補足: Falcon 仕様書 Table 3.2 の SamplerZ テストベクタ 16 本は、
   randombytes を**各チャンクごとに反転して**食わせないと再現しません
   （参照実装のテストハーネス `test.py` の `KAT_randbytes` が
   `bytes.fromhex(oc)[::-1]` としているのと同じ規約）。
   このファイルも同じ規約かどうかを知りたい。

## 出力の形式

- 一次資料（実際のファイルの中身）を最優先。
- `parameters.py` の該当箇所は**そのまま引用**してほしい。
- 見つからない場合は、round-1 / round-2 の提出パッケージや、
  Thomas Prest の Python 実装 <https://github.com/tprest/falcon.py> に
  同等のスクリプトがないかも調べてほしい。
```

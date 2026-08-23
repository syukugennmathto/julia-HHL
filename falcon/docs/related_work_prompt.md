# Deep Research 用プロンプト ― 関連研究の位置づけ

`docs/paper_outline.md` の §「投稿までにやること」の 1。
**これが終わるまで骨組みは仮**である。既知の結果を再発見している節が
ありうるからで、それを先に潰す必要がある。

文脈を知らないチャットに貼っても成立するよう自己完結させてある。

---

## プロンプト C ― FALCON の実装・再現性・サイドチャネルに関する先行研究

```
FALCON（NIST PQC、FN-DSA）の**実装**に関する先行研究を調べてください。
私は独立実装を書いて、参照実装の署名をビット単位で再現できることを示し、
その過程で仕様書が固定していない箇所をいくつか見つけました。
**それらが既に知られているかどうか**を知りたいのが目的です。

## 私が主張しようとしていること（この 4 つが新しいかを判定してほしい）

1. **仕様書だけを読んで書いた独立実装が、参照実装の署名を
   ビット単位で再現できる。** そして再現に必要だったのは
   浮動小数点の書き方を合わせることではなく、
   参照実装のサンプラ用 PRNG（SHAKE から導出される ChaCha20）の
   状態を外から固定する手段だった。

2. **仕様書と C 参照実装が「代数的に同値で丸めが異なる」箇所が 3 つある。**
   (a) 複素除算。C は `m = 1/(b_re²+b_im²)` を作って掛ける。
       仕様書に式はなく、Python 参照実装と Julia は言語の複素除算
       （Smith のアルゴリズム）を使う。
   (b) `LDL*` の `D11`。仕様書 Algorithm 8 と Python は
       `G11 − L10·adj(L10)·G00`（乗算 3 回）、
       C は `G11 − μ·adj(G01)`（μ = G01/G00、乗算 1 回）。
   (c) `ffSampling` の最下 2 段。C は logn==2 を手展開し、
       split/merge の捻り係数を `1/√2`, `1/√8` に畳んでいる。
   **そして測ると、この 3 つは署名に伝播する。ただし稀に。**
   n=512、鍵・メッセージ・PRNG 状態を揃えた 10 万署名のうち 5 本が
   食い違う（率 5e-5）。食い違うときは 512 係数のうち約 470 が違う。
   5 本すべてが整数近傍の中心 `mu`（例 221.00000000000006 対
   220.99999999999997）で、`SamplerZ` 冒頭の `floor(mu)` をまたいでいる。
   5 本すべてが 1024 回中 1023 回目か 1024 回目の呼び出し
   ＝位置 2n−2, 2n−1。
   一方 4 の定数（サンプラの**幅**を動かす）は、木の葉 512 枚中 444 枚が
   違うのに、同じ 10 万署名で差 0（95% 上限 3e-5）。
   境界は「差が大きいか小さいか」ではなく
   **差がどの量に届くか**（中心は `floor` を通り、幅は通らない）。

3. **仕様書 Table 3.3 の σ と σ_min は別の ε から導かれている。**
   σ は式 (2.13) を `ε ≤ 1/√(Q_s·λ)`（Q_s = 2^64、λ = 128/256）で
   評価した値（相対誤差 2e-16 で再現）。
   σ_min は `ε = 1/√(2^64·n³)` での Z の平滑化パラメータ
   （3e-13 で再現）。仕様書は後者の ε を印字していない。

4. **C 参照実装の `fpr_inv_sigma[logn]` は Table 3.3 の σ の
   正しく丸められた逆数ではない**（logn=9 で 1 ulp、logn=10 で 2 ulp）。
   参照実装は木の葉にこの表を掛けるので、`1/σ` を計算しても再現しない。

## 調べてほしい文献（最低限これらは読んで位置づけてほしい）

### 実装・再現性
- Thomas Pornin, "New Efficient, Constant-Time Implementations of Falcon"
  (IACR ePrint 2019/893)
- Pornin & Prest, "More Efficient Algorithms for the NTRU Key Generation
  using the Field Norm" (ePrint 2019/015)
- Pornin, "Improved (Again) Key Pair Generation for Falcon, BAT and Hawk"
  (ePrint 2025/1239)
- `pornin/c-fn-dsa` / `pornin/rust-fn-dsa`（FIPS 206 追随実装）
- PQClean / liboqs の Falcon 実装と、そこでの浮動小数点の扱い
- FALCON を固定小数点や整数で実装した研究
  （組込み向け、pqm4、Cortex-M4、FPGA、"fixed-point Falcon" 等）

### サンプラの等時間性
- Howe, Prest, Ricosset, Rossi, "Isochronous Gaussian Sampling:
  From Inception to Implementation" (PQCrypto 2020)
- Karmakar らの CDT / 表引きサンプラの定数時間実装

### 浮動小数点と側面攻撃
- **"Do Not Disturb a Sleeping Falcon" (ePrint 2024/1709)** ―
  浮動小数点誤差感度に基づく攻撃。**私の主張 2 と最も近い可能性がある。**
  この論文が「実装間の丸めの差」をどう扱っているかを詳しく。
- Fouque, Kirchner, Tibouchi, Wallet, Yu,
  "Key Recovery from Gram-Schmidt Norm Leakage in Hash-and-Sign
  Signatures over NTRU Lattices" (EUROCRYPT 2020, ePrint 2019/1180)
- Karabulut & Aysu, "Falcon Down: Breaking Falcon Post-Quantum Signature
  Scheme through Side-Channel Attacks" (DAC 2021)
- Guerreau, Martinelli, Ricosset, Rossi, "The Hidden Parallelepiped Is
  Back Again: Power Analysis Attacks on Falcon" (TCHES 2022)
- Zhang, Lin, Yu, Wang らの Falcon への電力・タイミング攻撃
- Fouque, Gajland, de Groote, Janneck, Kiltz,
  "A Closer Look at Falcon" (ePrint 2024/1769)

### 標準化
- NIST FIPS 206 の状態（2026 年 8 月時点で未発行）と、
  そこで予告されている変更（μ によるドメイン分離、randomized のみ、
  無限ノルム 840、公開鍵の NTT 化、**浮動小数点の演算順序の規定と
  FMA 禁止、署名の KAT との bit 一致要求**）
- 他の PQC 標準（ML-DSA/FIPS 204、SLH-DSA/FIPS 205）で
  「実装間のビット一致」がどう規定されているか

### 方法論
- Reparaz, Balasch, Verbauwhede, "Dude, is my code constant time?"
  (DATE 2017) と、その後の批判・改良
  （統計的検出力、雑音床、false positive の扱い）
- 暗号実装の性能比較における方法論（SUPERCOP、eBACS の測定手続き）

## 具体的に答えてほしいこと

1. **主張 1〜4 のそれぞれについて、既知か未知か。**
   既知なら誰がどこで述べているか（著者・論文・節）。
   「近いが違う」なら、どこがどう違うか。
2. 「**仕様書と参照実装の丸めの差**」を主題にした論文が既にあるか。
   FALCON 以外（ML-DSA, SLH-DSA, あるいは AES/SHA の時代）でも
   同種の研究があれば挙げてほしい。
3. **FALCON の独立実装で、参照実装と署名がビット一致すると
   主張しているものが既にあるか。** あればどの実装か。
4. 私の結果 2（3 箇所の書き方の差が 5e-5 で署名に伝播し、
   機構は `floor(mu)` の跨ぎ）は、既存のどの主張と整合し、
   どれと衝突するか。とくに **ePrint 2024/1709** について:
   (a) 摂動源として「ビルド構成の違い」以外を扱っているか。
       **代数的に同値な式の書き換え**を摂動源として挙げているか。
   (b) 報告している発散率の正確な値と、その測定条件。
   (c) 「幅を動かしても伝播しない」という対照実験があるか。
   この 3 点が、私の結果が再発見か否かを決める。
5. **投稿先の候補。** この内容（実装・再現性・標準化への含意が主、
   側面攻撃は付録）を受け入れる会議・論文誌はどこか。
   TCHES / CHES、PQCrypto、SAC、ACNS、JCEN、
   NIST PQC Standardization Conference のどれが適切か、
   それぞれの直近の採録傾向から判断してほしい。
6. **私の結果が「既知の再発見」で終わる危険が最も高いのはどれか。**
   忌憚なく。

## 出力の形式

- 一次資料（論文 PDF、実装のソース）を優先。ブログは導線としてのみ。
- 各主張に著者・年・会議/ePrint 番号・可能なら節番号。
- **「見つからなかった」も明記**してほしい。曖昧にぼかさないこと。
- 最後に、**「この 4 つの主張のうち、論文の主結果として立つのはどれか」**
  という判断を、理由付きで書いてください。
```

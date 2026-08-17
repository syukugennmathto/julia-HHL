# C 参照実装を dylib にする（macOS / Apple Silicon）

Dilithium のときと同じ「C 参照実装を `dylib` にして `ccall` で中間値を突き合わせる」
方法を FALCON でも使う。ただし **FALCON の参照実装は Dilithium ほど関数がフラットに
公開されていない**ので、その差分をここに書いておく。

> このドキュメントの記述は、実際に読んだファイルに基づく。
> ソースは <https://github.com/algorand/falcon>（round-3 参照実装のミラー。
> `Makefile` の著者表記は Thomas Pornin、ライセンスは MIT）。
> 行番号は 2026-08 時点の `master` のもの。

---

## 1. 取得とビルド

```sh
git clone https://github.com/algorand/falcon.git falcon-cref
cd falcon-cref
```

同梱の `Makefile` は `tests/` の実行ファイルを作るためのもので、共有ライブラリは
作らない（`Makefile` 冒頭、`OBJ = codec.o common.o deterministic.o falcon.o fft.o
fpr.o keygen.o rng.o shake.o sign.o vrfy.o`）。`ccall` 用にはこれを丸ごと 1 本の
`dylib` にする:

```sh
clang -O2 -fPIC -shared \
    -o libfalcon.dylib \
    codec.c common.c deterministic.c falcon.c fft.c fpr.c \
    keygen.c rng.c shake.c sign.c vrfy.c
```

- `-lm` は不要（`Makefile` のコメントいわく「native FPU を使う場合、x86 では
  通常不要」。arm64 も同様に libm への明示リンクなしで通る。リンクエラーが出たら
  足す）。
- `-O2` 以上が推奨（`Makefile` の `CFLAGS` は `-O3`）。ただし**突き合わせ用途では
  `-O0 -g` でビルドしたものも別名で用意しておくと lldb で中間値を覗ける**。

Universal binary が要るなら:

```sh
clang -O2 -fPIC -shared -arch arm64 -arch x86_64 -o libfalcon.dylib *.c
```

（`tests/` 以下の `main` を含むファイルを混ぜないよう `*.c` は使わず列挙するのが安全。）

---

## 2. Dilithium との最大の違い: `Zf()` プレフィックス

FALCON の内部関数は `inner.h` でこう定義されている:

```c
/* inner.h:290-295 */
#ifndef FALCON_PREFIX
#define FALCON_PREFIX   falcon_inner
#endif
#define Zf(name)             Zf_(FALCON_PREFIX, name)
#define Zf_(prefix, name)    Zf__(prefix, name)
#define Zf__(prefix, name)   prefix ## _ ## name
```

つまりソース中の `Zf(FFT)` は、実際のシンボル名としては **`falcon_inner_FFT`** に
なる。`ccall` で叩くときはこの展開後の名前を使う:

```julia
const LIBFALCON = "/path/to/libfalcon.dylib"

# void Zf(FFT)(fpr *f, unsigned logn);            inner.h:892
fft_cref!(f::Vector{Float64}, logn::Integer) =
    ccall((:falcon_inner_FFT, LIBFALCON), Cvoid, (Ptr{Float64}, Cuint), f, logn)
```

シンボルが見つからないときは、まず展開名を確認する:

```sh
nm -gU libfalcon.dylib | grep falcon_inner_ | head
```

`FALCON_PREFIX` を変えてビルドすれば別名にもできる。複数バージョンを同時に
ロードして差分を取りたくなったとき（round-3 と FIPS 206 draft の比較など）に効く。

---

## 3. 落とし穴 (1): 既定では浮動小数点が**エミュレーション**される

> **訂正**: この節の初版には「Apple Silicon では `FALCON_FPNATIVE` が選ばれる」
> と書いていた。`inner.h:189-201` の既定ロジックだけを見ればそうなるが、
> **`config.h:87` が `#define FALCON_FPEMU 1` で無条件に上書きしている**ので、
> 素のままビルドすると**エミュレーション版になる**。実際にビルドして確認済み
> （`docs/debug_log.md` #012）。`-DFALCON_FPNATIVE=1` を付けても効かない
> （`config.h` が後勝ち）。

`config.h:64-86` は、その理由を *** CRITICAL SECURITY WARNING *** 付きで
書いている。要約:

- ネイティブ FPU やコード最適化はわずかな差異を生み、**決定的署名**を壊しうる
- 署名の非決定性は**壊滅的なセキュリティ障害**につながる。同一メッセージに
  対する異なる 2 つの署名から、任意メッセージの偽造が可能になりうる
- したがって `FALCON_FPEMU` を有効のままにすることを**強く推奨**する

つまり FALCON では「Float64 の結果がプラットフォーム間で bit 単位に
再現すること」自体がセキュリティ性質である。詳細は `docs/math/05_fft.md`。

### 2 つのモードでの `fpr`

| モード | `fpr` の型 | 定義 |
|:---|:---|:---|
| `FALCON_FPEMU=1`（既定） | `uint64_t` | `fpr.h:108` |
| `FALCON_FPNATIVE=1` | `struct { double v; }` | `fpr.h:510` |

**エミュレーション版の `uint64_t` は IEEE-754 binary64 の bit パターン
そのもの**である（`fpr_gm_tab` の値をデコードして確認済み）。
だから `ccall` する側では、どちらのモードでも 8 バイトの塊として扱い、
Julia 側では `Float64` として reinterpret すれば通る。

- **配列として渡す分には両モードとも問題ない。** サイズ 8・アライメント 8 で
  一致するので、`Ptr{Float64}` で `fpr *` を受ける関数は正しく動く。
- **スカラを値渡しする関数は注意。** FPNATIVE の
  `struct { double v; }` は arm64 では HFA として浮動小数点レジスタで渡されるが、
  規約に寄りかかった書き方になる。値渡し API が要るなら C 側に薄いラッパを:

```c
/* shim.c -- 突き合わせ専用 */
#include "inner.h"
double shim_fpr_add(double a, double b) {
    fpr x, y, z;
#if FALCON_FPEMU
    memcpy(&x, &a, 8); memcpy(&y, &b, 8);
#else
    x.v = a; y.v = b;
#endif
    z = fpr_add(x, y);
    double d;
#if FALCON_FPEMU
    memcpy(&d, &z, 8);
#else
    d = z.v;
#endif
    return d;
}
```

## 3.5 落とし穴 (2): `fpr_double()` は「double へ変換」ではない

```c
/* fpr.h:682 (FPNATIVE) */
static inline fpr fpr_double(fpr x) { return FPR(x.v + x.v); }
```

**2 倍する関数**である。`fpr_half` の対になっている。
`double` への変換だと思って使うと、ダンプが全部ゼロになるか、
もっと悪いことに 2 倍された値が返る。実際に踏んだ
（`docs/debug_log.md` #011）。

## 3.6 落とし穴 (3): FFT の**表現**が違う

`Zf(FFT)` の出力は、我々の（および Python 参照実装の）表現とは違う。
`inner.h:875-884` にこうある:

> 実多項式は N 個の `fpr` の配列で表現される。実多項式の FFT 表現は
> **N/2 個の複素要素**を含み、各々は実部・虚部の 2 つの実数として格納される。

つまり:

1. **半分しか持たない。** 実多項式の FFT は共役対称なので後半は冗長。
2. **実部と虚部を分けて置く。** `n/2` 個の実部の後に `n/2` 個の虚部。
3. **順序が違う。** 前半 `n/2` 個の中での並びが、我々の並びと
   **グレイコード置換**でずれている。C の添字 `k` は我々の添字
   `k XOR (k >> 1)`（0 始まり）に対応する。

3 番目は `inner.h` に書かれていない（"see falcon-fft.c for details" とあるだけ）。
実際にビルドして n = 4〜1024 で突き合わせて determine した。

Julia 側には `Falcon.from_c_fft` / `Falcon.to_c_fft` として実装済み。
また、C から取った実際のベクタが `test/vectors/fft_c_kat.jl` にコミットして
あるので、**dylib をビルドしなくてもモジュール 5 は C と突き合わせ済み**である。
再生成したい場合のドライバは `scripts/cref_fft_dump.c`。

## 4. 突き合わせに使える内部関数（`inner.h` の行番号つき）

`Zf(...)` は `falcon_inner_...` に読み替えること。

| モジュール | C 関数 | 宣言 |
|:---|:---|:---|
| shake | `Zf(i_shake256_init/inject/flip/extract)` | inner.h:428-434 |
| PRNG | `Zf(prng_init)`, `Zf(prng_refill)`, `Zf(prng_get_bytes)` | inner.h:807-818 |
| hash | `Zf(hash_to_point_vartime)`, `Zf(hash_to_point_ct)` | inner.h:536, 547 |
| verify | `Zf(is_short)`, `Zf(is_short_half)` | inner.h:556, 568 |
| verify | `Zf(verify_raw)`, `Zf(verify_recover)` | inner.h:592, 659 |
| NTT | `Zf(to_ntt_monty)` | inner.h:579 |
| 鍵 | `Zf(compute_public)`, `Zf(complete_private)`, `Zf(is_invertible)` | inner.h:604, 618, 628 |
| FFT | `Zf(FFT)`, `Zf(iFFT)` | inner.h:892, 902 |
| FFT 領域演算 | `Zf(poly_add/sub/neg/adj_fft/mul_fft/muladj_fft/mulselfadj_fft)` | inner.h:908-944 |

**ただし「関数が一覧にある」ことと「そのまま比較できる」ことは別**である
（3.6 節）。プロトタイプは表現を教えてくれない。実際に突き合わせる前に、
必ず小さい n で 1 回ダンプして表現を確認すること。

モジュール 3（`poly.jl`）に対応する C 関数が無いのは、参照実装が素朴な畳み込みを
持っていないから（常に NTT か FFT を使う）。だから**モジュール 3 の突き合わせ相手は
Python 参照実装**にしてある（`scripts/pyref/ntt.py` の `mul_zq`）。
モジュール 4（`ntt.jl`）も同様で、C の `Zf(to_ntt_monty)` はビット反転順かつ
Montgomery 形式なので Python 側を正典にした（`docs/debug_log.md` #007）。

---

## 5. KAT の入手

`deterministic.c` が固定シードからの決定的な鍵生成・署名を提供している
（`Makefile` の `OBJ` に含まれる）。KAT を作るときはここを入口にする。
Python 側の KAT は `tprest/falcon.py` の `scripts/sign_KAT.py` などにある
（未取得。モジュール 10 で必要になったら取る）。

---

## 6. `ccall` の作法（Dilithium のときと同じだが再掲）

- 配列は `Ptr{T}` で渡す。Julia の `Vector` は連続メモリなので `pointer` 不要、
  そのまま渡してよい。GC には `ccall` が引数を保持するので追加の `GC.@preserve` は
  基本不要。
- **C 側が破壊的に書き換える関数が多い**（`Zf(FFT)` は in-place）。突き合わせでは
  必ず入力のコピーを渡すこと。ここを忘れると「2 回目から値が違う」という
  再現性のないバグに化ける。
- `unsigned logn` は `Cuint`。`n` ではなく `log2(n)` を渡す。**これが FALCON で
  一番よくやる引数ミス**で、`logn=512` を渡すとシフト量が壊れて静かに segfault
  するか、もっと悪いことに壊れた値を返す。

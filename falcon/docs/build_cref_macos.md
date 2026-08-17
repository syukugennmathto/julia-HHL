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

## 3. 落とし穴: `fpr` は `double` そのものではない

Apple Silicon（`__ARM_FP` が定義される）では、`inner.h:189-201` の分岐により
`FALCON_FPNATIVE = 1` が選ばれる。このとき `fpr.h:510`:

```c
typedef struct { double v; } fpr;
```

**`double` を包んだ struct** であって `double` ではない。実務上の帰結:

- **配列として渡す分には問題ない。** `struct { double v; }` のメモリ表現は
  `double` と同一（サイズ 8、アライメント 8）なので、`Ptr{Float64}` で
  `fpr *` を受ける関数は正しく動く。上の `Zf(FFT)` がこれ。
- **スカラを値渡しする関数は注意。** arm64 の呼び出し規約では「double 1 個だけの
  struct」は HFA として浮動小数点レジスタで渡されるので結果的に `Cdouble` と
  一致するが、規約に寄りかかった書き方になる。値渡し API を叩く必要が出たら、
  C 側に薄いラッパを 1 本足して `double` で受け渡しするほうが安全:

```c
/* shim.c -- 突き合わせ専用の薄いラッパ */
#include "inner.h"
double shim_fpr_add(double a, double b) {
    fpr x, y, z;
    x.v = a; y.v = b;
    z = fpr_add(x, y);
    return z.v;
}
```

- `FALCON_FPEMU=1` でビルドすると `fpr` は `uint64_t`（`fpr.h:108`）になり、
  浮動小数点演算が整数演算でエミュレートされる。**これは定数時間実装のための
  ものであり、原稿の主題そのもの**なので、後で `FALCON_FPEMU=1` 版も別名で
  ビルドして「native FPU 版と bit 単位で一致するか」を見ておく価値がある
  （仕様が Float64 の丸めまで固定していることの実地確認になる）。

---

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

モジュール 5（`fft.jl`）以降はこの表の関数で 1 対 1 に突き合わせられる。
モジュール 3（`poly.jl`）に対応する C 関数が無いのは、参照実装が素朴な畳み込みを
持っていないから（常に NTT か FFT を使う）。だから**モジュール 3 の突き合わせ相手は
Python 参照実装**にしてある（`scripts/pyref/ntt.py` の `mul_zq`）。

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

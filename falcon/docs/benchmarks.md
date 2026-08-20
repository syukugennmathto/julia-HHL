# C 参照実装 vs. この Julia 実装 ― 速度対決

> 生データと再現手順のみ。設計上の解釈は `docs/debug_log.md` #031 / #032 に書く。

## 測定条件

| | |
|:---|:---|
| CPU | x86-64（コンテナ内、`fma` あり）。**Apple Silicon ではない** |
| Julia | 1.12.7（`docs/debug_log.md` #013 の手順で入手） |
| C | round-3 参照実装（algorand/falcon ミラー）、`gcc -O2` |
| 単位 | 1 回あたりのミリ秒、**中央値** |
| 反復 | keygen 20 / sign 200 / verify 200（C）、keygen 5 / sign 30 / verify 100（Julia） |

再現:

```sh
# C 側（cbuild/ に参照実装の .c を置いた状態で）
cc -O2 -o cref_bench cref_bench.c \
   codec.c common.c falcon.c fft.c fpr.c keygen.c rng.c shake.c sign.c vrfy.c -lm
cc -O2 -DFALCON_FPEMU=0 -DFALCON_FPNATIVE=1 -o cref_bench_nat cref_bench.c ...同上...
./cref_bench     9 20 200 200 > c_emu_512.txt
./cref_bench_nat 9 20 200 200 > c_nat_512.txt

# Julia 側
julia --project=falcon falcon/scripts/bench.jl 9 5 30 100 > jl_512.txt

# 突き合わせ
julia falcon/scripts/bench_compare.jl c_emu_512.txt c_nat_512.txt jl_512.txt
```

**中央値を採る理由**: 鍵生成は `ntru_solve` が解けるまで再試行するので
分布に長い右裾がある。平均は代表値にならない
（C の実測でも min 10.9 ms / median 13.0 ms / mean 16.6 ms と開く）。

**Julia 側は計測前に全操作を 1 回空回ししている**。最初の呼び出しは
実行ではなくコンパイルであり、混ぜると 100 倍単位で嘘の数字になる。
GC は止めていない ― 割り当てはこの実装の実コストであり、隠すのは不誠実。

**比較の単位を合わせること**: C の `falcon_sign_dyn` は毎回 Falcon 木を
作り直し、`falcon_sign_tree` は展開済み鍵を使う。我々の `falcon_sign` は
展開済み鍵を取るので、**対応するのは `sign_tree`**。
`sign_dyn` と比べると木の構築コスト分だけ我々に下駄を履かせることになる。

## FALCON-512

| op | C (FPEMU) | C (native FP) | Julia | Julia / C(FPEMU) |
|:---|---:|---:|---:|---:|
| `keygen` | 13.0 | 6.35 | **9610** | 738x |
| `expand_privkey` | 1.77 | 0.081 | **0.980** | **0.6x** |
| `sign_dyn` | 4.43 | 0.311 | （該当なし） | — |
| `sign_tree` | 2.04 | 0.194 | **2.05** | **1.0x** |
| `verify` | 0.030 | 0.029 | **0.458** | 15.3x |

### 読み方

**1. 署名は C の FPEMU 版と完全に互角**（2.05 ms 対 2.04 ms）。
鍵展開に至っては Julia のほうが速い（0.98 ms 対 1.77 ms）。

これは Julia が速いというより、**FPEMU が高い**という話である。
C 側はネイティブ FPU を使えば署名が 10 倍速くなる（2.04 → 0.194 ms）。
つまり参照実装は、移植性・再現性のために**実行速度の 1 桁を支払っている**。
我々はハードウェアの double をそのまま使っているので、
その 1 桁を払っていない ― 代わりに、払わなかったぶんだけ
「どのマシンでも同じ署名が出る」という保証を失っている
（`docs/debug_log.md` #031）。

**この 2 つは同じコインの裏表であり、性能表だけを見ても意味が取れない。**

**2. 検証が 15 倍遅いのは浮動小数点のせいではない。**
`falcon_verify` は**整数演算だけ**である（`src/falcon.jl` の docstring）。
だから「Julia が遅いのは FP まわりだろう」は仮説として外れている（#032）。

内訳を測った（n=512、中央値 ms）:

| 段階 | 時間 | 割合 |
|:---|---:|---:|
| `decode_signature` | 0.0065 | 2% |
| `hash_to_point` | 0.026 | 7% |
| `polymulq`（学校算法） | 0.126 | 33% |
| **`sqnorm`（`BigInt`）** | **0.186** | **49%** |
| 合計 | 0.345 | |
| `falcon_verify` 全体 | 0.380 | |

**一番高いのは多項式乗算ではなく `sqnorm` である。**
`BigInt` の割り当てが 1024 係数ぶん走る。C は `int32` の
積和を回すだけなので、ここは丸ごと余分なコストになっている。

そして、参考のために測った `polymulq_ntt` が

| | 時間 |
|:---|---:|
| `polymulq`（学校算法, O(n²)） | 0.126 ms |
| `polymulq_ntt`（NTT 経由, O(n log n)） | **0.938 ms** |

**NTT のほうが 7.4 倍遅い。**
「`verify` を NTT に差し替えれば速くなる」は完全に間違いだった。
n=512 では O(n²) の 262144 回の積和のほうが、
再帰的な `merge_ntt`/`split_ntt` が各段で配列を割り当てるコストより安い。
漸近計算量が勝つ次数にまだ届いていない
（そして FALCON には n=1024 より上が無い）。

**3. 鍵生成の 738 倍が本丸。**
`ntru_solve` の再帰下降が `BigInt` を使っている。
C は 31 bit リムの RNS（剰余系）で同じ計算をやり、割り当てもしない。
これは定数倍ではなく**アルゴリズムの差**であり、
他の 3 項目（0.6x / 1.0x / 15x）と桁が違うことがその証拠になっている。

なお `BigInt` を使っているのは怠慢ではない。
n=512 の降下の底では係数が数千 bit になり、`Int128` でも足りない
（`docs/math/06_ntrugen.md`）。速くするなら「BigInt をやめる」のではなく
「C と同じ RNS を実装する」ことになる。

## FALCON-1024

| op | C (FPEMU) | C (native FP) | Julia | Julia / C(FPEMU) |
|:---|---:|---:|---:|---:|
| `keygen` | 55.7 | 21.1 | **88007** | 1580x |
| `expand_privkey` | 3.94 | 0.175 | **1.81** | **0.5x** |
| `sign_dyn` | 9.93 | 0.648 | （該当なし） | — |
| `sign_tree` | 4.29 | 0.489 | **4.12** | **1.0x** |
| `verify` | 0.061 | 0.064 | **1.06** | 17.5x |

**FALCON-1024 は keygen / sign / verify がすべて通る**
（`test/test_falcon1024.jl`。この表もその実行そのものである）。

n=512 との比を見ると、性質が 512 の表と揃っている:

| | Julia の 1024/512 比 | C(FPEMU) の 1024/512 比 |
|:---|---:|---:|
| `keygen` | 9.2x | 4.3x |
| `sign_tree` | 2.0x | 2.1x |
| `verify` | 2.3x | 2.0x |

署名と検証はどちらの実装でもほぼ 2 倍 ―
n が 2 倍で、支配項が係数ごとの定数作業だからである。
**鍵生成だけが 9.2 倍**（C は 4.3 倍）で、
Julia のほうが次数の増加に対して**さらに悪化する**。
これは `ntru_solve` の降下が深くなり、
底の係数の bit 長も伸びて `BigInt` の負担が二重に増えるためで、
「定数倍ではなくアルゴリズムの差」という #032 の読みと整合している。

## この表を回帰テストにしてはいけない

しきい値を主張するテストは書いていない。ベンチマークは実行環境の負荷で
動くので、CI に入れれば必ず偽陽性を出す。
性能の変化を追いたいなら、この文書の数値を**手で**更新して
差分をコミットに残すのが正しい。

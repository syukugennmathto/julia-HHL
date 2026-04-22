# julia-HHL

「線形方程式の解き方を、量子で」のコードスニペットを **Python (Qiskit)** と
**Julia (Yao.jl)** の両方で実装・実行し、動作確認とベンチマーク比較を行った記録。

## 対象問題

2×2 エルミート線形方程式:

```
A = [[1, -1/3], [-1/3, 1]],  b = [1, 0]
古典解: x = [1.125, 0.375]  (正規化後 [0.949, 0.316])
```

## ディレクトリ

```
python/     Qiskit 実装
  ch2_basics.py       基礎ゲート (H, X, Z, 測定)
  ch3_encoding.py     基底/テンソル積/量子カーネル/QPE 骨格
  ch4_hhl.py          HHL 2x2 本体 + ノイズシミュレーション
  bench_hhl.py        細粒度ベンチマーク

julia/      Yao.jl 実装
  ch5_hhl.jl          HHL 2x2 本体
  bench_hhl.jl        細粒度ベンチマーク

results/    実行結果
  comparison.md       総合比較レポート
  *.json              各ベンチマークの数値
```

## 実行方法

### Python

```bash
python3 -m venv .venv
.venv/bin/pip install qiskit qiskit-aer
.venv/bin/python python/ch4_hhl.py
.venv/bin/python python/bench_hhl.py
```

### Julia

```bash
julia -e 'using Pkg; Pkg.activate("julia_env"; shared=false); Pkg.add(["Yao","BenchmarkTools","JSON"])'
julia julia/ch5_hhl.jl
julia julia/bench_hhl.jl
```

## 主な結果

- HHL 2x2 で Python/Julia 両実装とも **fidelity 0.900** で一致 (`n=2` 精度制約)。
- 純シミュ計算 (run-only) では Julia/Yao.jl が Python/Qiskit より **20〜40倍** 高速。
- end-to-end (build + transpile + run) では Python 側の transpile オーバーヘッド
  (~110ms) が支配的で、**約 6000 倍** の差。

詳細は [`results/comparison.md`](results/comparison.md)。

## 環境

- macOS (arm64)
- Python 3.9.6 / Qiskit 2.2.3 / Qiskit-Aer 0.17.2
- Julia 1.12.4 / Yao.jl 0.9.x

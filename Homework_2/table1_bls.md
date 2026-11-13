# Table 1 (Problem 5, Task 3) — draft

Source: `table1` from `q5.R`, R 4.3.3, run 2026-10-21 in WSL Ubuntu.
Runtimes vary a little between runs, so rerun before submitting if the report must match a specific run.

**Table 1.** Cost of GD with backtracking line search for every distinct (c, β) setting used in Problem 5. All runs start at x₀ = x_init with trial step γ = 1 reset at every iteration. The shared setting (c, β) = (0.1, 0.5) appears only once. The objective-evaluation count includes L(x₀) and every Armijo trial, accepted or rejected. Accepted trial values are cached and reused as L(xₙ). Runtime is the median wall-clock time of 21 timed runs after one untimed warm-up, excluding CSV loading, the precomputation of G, v and C, and the SVD benchmark. IQR is the interquartile range of those 21 timings.

| c    | β   | L(x₃₀) − L_opt  | GD updates | Objective evals | Backtracking reductions | Median runtime (ms) | IQR (ms) |
|------|-----|-----------------|-----------:|----------------:|------------------------:|--------------------:|---------:|
| 0.1  | 0.1 | 1.352 × 10⁻¹    | 30 | 116  | 85   | 1.16  | 0.30 |
| 0.1  | 0.2 | 5.157 × 10⁻³    | 30 | 138  | 107  | 1.07  | 0.06 |
| 0.1  | 0.5 | 1.132 × 10⁻⁴    | 30 | 267  | 236  | 1.89  | 0.13 |
| 0.1  | 0.9 | 6.285 × 10⁻²    | 30 | 1658 | 1627 | 10.82 | 0.79 |
| 0.01 | 0.5 | 6.237 × 10⁻⁴    | 30 | 269  | 238  | 1.92  | 0.09 |

Check: evals = 1 (for L(x₀)) + 30 (accepted trials) + reductions, which holds for every row.

## LaTeX version

```latex
\begin{table}[h]
\centering
\caption{Cost of GD with backtracking line search for every distinct $(c,\beta)$
setting used in Problem~5. The shared setting $(c,\beta)=(0.1,0.5)$ appears once.
Objective evaluations include $L(x_0)$ and every Armijo trial, accepted or rejected;
accepted trial values are cached as $L(x_n)$. Runtime is the median wall-clock time
of 21 timed runs after one untimed warm-up, excluding CSV loading, precomputation of
$G$, $v$, $C$, and the SVD benchmark. IQR is the interquartile range of the 21 timings.}
\label{tab:bls}
\begin{tabular}{cccrrrrr}
\hline
$c$ & $\beta$ & $L(x_{30})-L_{\mathrm{opt}}$ & Updates & Obj.\ evals & Reductions & Median (ms) & IQR (ms) \\
\hline
0.1  & 0.1 & $1.352\times10^{-1}$ & 30 & 116  & 85   & 1.16  & 0.30 \\
0.1  & 0.2 & $5.157\times10^{-3}$ & 30 & 138  & 107  & 1.07  & 0.06 \\
0.1  & 0.5 & $1.132\times10^{-4}$ & 30 & 267  & 236  & 1.89  & 0.13 \\
0.1  & 0.9 & $6.285\times10^{-2}$ & 30 & 1658 & 1627 & 10.82 & 0.79 \\
0.01 & 0.5 & $6.237\times10^{-4}$ & 30 & 269  & 238  & 1.92  & 0.09 \\
\hline
\end{tabular}
\end{table}
```

# Table 2 (Problem 6, Task 2) — draft

Source: `table2` from `q6.R`, R 4.3.3, WSL Ubuntu. Final gaps, update counts and evaluation counts are deterministic. Runtimes come from the timings saved in `timings_problem6.csv`. A rerun on 2026-10-23 gave the same medians to within about 10%, but one BLS run had outliers (IQR ≈ 2.4 ms), so the committed timings were kept. Rerun before submitting if the report must match a specific run.

**Table 2.** Cost of the three GD configurations compared in Problem 6. All runs start at x₀ = x_init and perform 30 updates (the gradient-norm tolerance 10⁻⁸ was never reached). Objective evaluations count only those the algorithm needs. FSS and ELS need none. For BLS the count includes L(x₀) and every Armijo trial, accepted or rejected, with accepted trial values cached as L(xₙ). Runtime is the median wall-clock time of 21 timed runs after one untimed warm-up. **Runtimes exclude CSV loading, the common precomputation of G = AᵀA, v = Aᵀy and C, and the SVD benchmark (L_opt, x_opt).** IQR is the interquartile range of the 21 timings.

| Method                 | L(x₃₀) − L_opt | GD updates | Objective evals | Median runtime (ms) | IQR (ms) |
|------------------------|----------------|-----------:|----------------:|--------------------:|---------:|
| FSS, γ = 1/C (p = 1)   | 9.808 × 10⁻¹   | 30 | 0   | 0.152 | 0.017 |
| ELS                    | 4.436 × 10⁻²   | 30 | 0   | 0.257 | 0.020 |
| BLS, c = 0.1, β = 0.5  | 1.132 × 10⁻⁴   | 30 | 267 | 1.893 | 0.046 |

Check: the BLS row matches the (0.1, 0.5) row of Table 1 (267 evals = 1 + 30 + 236 reductions).

## Work not captured by the evaluation count (for the text under Table 2)

Per iteration, with N = 100 and K = 50:

| Method | Gradient g = 2(Gx − v) | Step-size work | Matrix–vector products per iteration |
|--------|------------------------|----------------|--------------------------------------|
| FSS | one K×K product | none (γ fixed) | 1 (K×K) |
| ELS | one K×K product | one extra K×K product Gg for γ = ‖g‖² / (2gᵀGg) | 2 (K×K) |
| BLS | one K×K product | one N×K product A x per Armijo trial | 1 (K×K) + (1 + backtracks) (N×K) |

Over 30 iterations that is 30 K×K products for FSS, 60 for ELS, and 30 K×K + 266 N×K for BLS. This ordering matches the runtimes. ELS has zero objective evaluations, yet it costs about 1.7× FSS because of the extra Gg product.

## LaTeX version

```latex
\begin{table}[h]
\centering
\caption{Summary of the computational cost of FSS, ELS, and BLS.}
\label{tab:comparison}
\begin{tabular}{lcrrrr}
\hline
Method & $L(x_{30})-L_{\mathrm{opt}}$ & Updates & Obj.\ evals & Median (ms) & IQR (ms) \\
\hline
FSS, $\gamma = 1/C$            & $9.808\times10^{-1}$ & 30 & 0   & 0.152 & 0.017 \\
ELS                            & $4.436\times10^{-2}$ & 30 & 0   & 0.257 & 0.020 \\
BLS, $c=0.1$, $\beta=0.5$      & $1.132\times10^{-4}$ & 30 & 267 & 1.893 & 0.046 \\
\hline
\end{tabular}
\end{table}
```

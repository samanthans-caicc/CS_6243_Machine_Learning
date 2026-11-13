# HW2 Key Discussion Points: Gradient Descent for Least Squares

Draft talking points for the report, organized by problem. Numbers come from running
`q0.R`–`q6.R` (R 4.3.3, WSL Ubuntu, 2026-10-22). Runtimes change slightly from run to run,
so take the final values from the run you actually report.

---

## Setup and conventions (Problem 0)

- Data dimensions: $A \in \mathbb{R}^{100 \times 50}$, $y \in \mathbb{R}^{100}$, $x_{\text{init}} \in \mathbb{R}^{50}$, so $N = 100$, $K = 50$.
- Objective and precomputed quantities:
  $$L(x) = \|Ax - y\|_2^2, \qquad G = A^\top A, \qquad v = A^\top y, \qquad g(x) = 2(Gx - v). \tag{1}$$
- GD update, with up to 30 updates starting from $x_0 = x_{\text{init}}$:
  $$x_n = x_{n-1} - \gamma_n\, g(x_{n-1}), \qquad n = 1, \dots, 30. \tag{2}$$
- Gradient-norm stopping tolerance: $\|g(x_{n-1})\|_2 \le 10^{-8}$, checked **before** each step. No method triggered it, so every run performed all 30 updates.
- The runs use no randomness, so no seed is needed. Arithmetic is IEEE double precision (R's default numeric type).
- **Flag:** the PDF says "Use Python" (Problems 2–5) and suggests `time.perf_counter()`. The code is in R and times with `Sys.time()`. Either port the code or state and justify the deviation in the report.

---

## Problem 1: Derivations

**Gradient.** Expand the quadratic:
$$L(x) = x^\top A^\top A x - 2y^\top A x + y^\top y \;\Longrightarrow\; g(x) = 2A^\top(Ax - y) = 2(Gx - v). \tag{3}$$

**Hessian.** The Hessian is constant (it does not depend on $x$):
$$H = \nabla^2 L(x) = 2A^\top A. \tag{4}$$

**Lipschitz constant.** Since $g(u) - g(v) = H(u - v)$,
$$\|g(u) - g(v)\|_2 \le \|H\|_2\,\|u - v\|_2, \qquad C = \|H\|_2 = \lambda_{\max}(H) = 2\sigma_1^2. \tag{5}$$
This $C$ is the *smallest* valid constant because equality holds when $u - v$ is the top eigenvector of $H$.

**Set of LS minimizers.** With the rank-reduced SVD $A = U_r \Sigma_r V_r^\top$, where $r = \operatorname{rank}(A)$:
$$\mathcal{X}^\star = \left\{ V_r \Sigma_r^{-1} U_r^\top y + z \;:\; z \in \mathcal{N}(A) = \operatorname{range}(V_\perp) \right\}, \tag{6}$$
$$L_{\text{opt}} = \|y\|_2^2 - \|U_r^\top y\|_2^2 = \|(I - U_r U_r^\top)\,y\|_2^2. \tag{7}$$

**Minimum-norm minimizer.** $x_{\text{opt}} = V_r \Sigma_r^{-1} U_r^\top y = A^{+} y$. Because $x_{\text{opt}} \in \operatorname{range}(V_r) \perp z$, Pythagoras gives
$$\|x\|_2^2 = \|x_{\text{opt}}\|_2^2 + \|z\|_2^2 \ge \|x_{\text{opt}}\|_2^2. \tag{8}$$
This argument also covers rank-deficient $A$. When $r = K$, $\mathcal{N}(A) = \{0\}$ and the minimizer is unique.

---

## Problem 2: Numerical results

| Quantity | Value |
|---|---|
| $C = \lambda_{\max}(H) = 2\sigma_1^2$ | $565.0373$ |
| $\lambda_{\min}^{+}(H) = 2\sigma_K^2$ | $25.5747$ |
| $\sigma_1,\ \sigma_K$ | $16.808,\ 3.576$ |
| Condition number $\kappa(H) = \lambda_{\max}/\lambda_{\min}$ | $\approx 22.09$ (that is, $4.70^2$) |
| Rank tolerance $\max(N,K)\,\varepsilon\,\sigma_1$ | $3.73 \times 10^{-13}$ |
| Numerical rank $r$ | $50$ (full column rank) |
| $L_{\text{opt}}$ (residual form and formula (7) agree) | $37.5113760783912$ |
| $\|x_{\text{opt}}\|_2$ | $0.9736$ |

- **Key point:** $A$ has full column rank, so $H \succ 0$, the minimizer is unique, and $x_{\text{opt}}$ is reachable by GD. The null-space discussion in Problem 3 is therefore theoretical for this dataset.

---

## Problem 3: Fixed step size $\gamma = p/C$

**Theory.** Write the error as $e_n = x_n - x_{\text{opt}}$. In the eigenbasis of $H$, each component evolves independently:
$$e_n = (I_K - \gamma H)\, e_{n-1} \;\Longrightarrow\; q_i^\top e_n = (1 - \gamma\lambda_i)^n\, q_i^\top e_0. \tag{9}$$
Convergence from **every** initialization requires $|1 - \gamma\lambda_i| < 1$ for all $\lambda_i > 0$:
$$0 < \gamma < \frac{2}{\lambda_{\max}} = \frac{2}{C} \quad\Longleftrightarrow\quad 0 < p < 2. \tag{10}$$
- The lecture's $\gamma \le 1/C$ ($p \le 1$) is a conservative *sufficient* condition. The exact LS range in (10) is twice as wide.
- If $H$ were singular, then for $\lambda_i = 0$ the factor is $(1 - 0) = 1$. The null-space component of $e_0$ never changes, so GD converges to the minimizer closest to $x_0$, not necessarily to $x_{\text{opt}}$. This does not occur here because $r = K$.

**Observed results.** The contraction factor is $\rho(p) = \max_i |1 - \gamma\lambda_i|$.

| $p$ | $\rho(p)$ | $L(x_{30}) - L_{\text{opt}}$ | $\|x_{30} - x_{\text{opt}}\|_2$ | Behavior |
|---|---|---|---|---|
| 0.1 | 0.9955 | $355.69$ | 2.691 | Too timid: monotone but very slow |
| 0.5 | 0.9774 | $10.90$ | 0.703 | Monotone, slow |
| 1.0 | 0.9547 | $0.981$ | 0.238 | Top component removed in **one** step ($1 - 1 = 0$) |
| 1.5 | 0.9321 | $0.134$ | 0.095 | **Best**: top component oscillates in sign but decays ($\times(-0.5)$ per step) |
| 2.0 | 1.0000 | $253.80$ | 0.949 | Boundary case: top component flips sign forever ($\times(-1)$ per step) |

- **Objective vs. iterates:**
  - For $p = 1.5$ the objective decreases monotonically, yet the iterate component along $q_1$ **alternates sign**: $-0.948 \to 0.474 \to -0.237 \to \cdots$. This oscillation cannot be seen in $L(x_n)$ because $L$ depends on $(q_1^\top e_n)^2$.
  - For $p = 2$, $L(x_n)$ plateaus at a high value while the iterates bounce between two points along $q_1$ ($\pm 0.948$). All other components still decay, which explains why the objective drops at first and then stalls.
- **Guarantee vs. this initialization:** $p = 2$ is outside the open interval in (10), so convergence is not guaranteed. Here it fails because $q_1^\top e_0 \ne 0$. A special $x_0$ with $q_1^\top e_0 = 0$ would still converge.
- The step that minimizes the worst-case contraction factor is
  $$\gamma^\star = \frac{2}{\lambda_{\max} + \lambda_{\min}}, \qquad \rho^\star = \frac{\kappa - 1}{\kappa + 1} \approx 0.913, \qquad p^\star = \frac{2C}{C + \lambda_{\min}} \approx 1.91. \tag{11}$$
  This explains why $p = 1.5$ beats $p = 1$.

---

## Problem 4: Exact line search

**Closed-form step.** Minimizing $\phi(\gamma) = L(x - \gamma g)$ over $\gamma$ gives
$$\gamma_n = \frac{\|g_{n-1}\|_2^2}{g_{n-1}^\top H\, g_{n-1}} = \frac{\|g\|_2^2}{2\,g^\top G g}. \tag{12}$$

**Bounds.** The Rayleigh quotient of $H$ lies in $[\lambda_{\min}, \lambda_{\max}]$, so
$$\frac{1}{C} \;\le\; \gamma_n \;\le\; \frac{1}{\lambda_{\min}^{+}}, \qquad\text{i.e.}\qquad 0.00177 \le \gamma_n \le 0.0391. \tag{13}$$

**Observations.**
- Observed $\gamma_n C$ ranges from $1.52$ to $1.90$, so every ELS step is **larger than $1/C$**. ELS adapts to the curvature *along $g$*, which is less than the worst-case curvature $C$.
- $\gamma_n$ fluctuates from step to step. This is the classic **zig-zag**: consecutive gradients are orthogonal ($g_n^\top g_{n-1} = 0$), and the step alternates between directions with different curvature.
- Final gap: $L(x_{30}) - L_{\text{opt}} = 0.0444$. This beats FSS with $p = 1$ ($0.981$), but ELS is only locally optimal per step. Its worst-case rate is
  $$L(x_n) - L_{\text{opt}} \le \left(\frac{\kappa - 1}{\kappa + 1}\right)^{2n} \big(L(x_0) - L_{\text{opt}}\big). \tag{14}$$
- **Safeguard:** if $g^\top G g \le 10^{-300}$ (no curvature along $g$), stop and report rather than divide by zero. This safeguard never triggered.

---

## Problem 5: Backtracking line search (Armijo)

**Acceptance rule.** For this quadratic, the Armijo condition can be solved exactly:
$$L(x - \gamma g) = L(x) - \gamma\|g\|^2 + \tfrac{\gamma^2}{2}\, g^\top H g \le L(x) - c\gamma\|g\|^2 \;\Longleftrightarrow\; \gamma \le \bar\gamma_c := \frac{2(1-c)\,\|g\|^2}{g^\top H g}. \tag{15}$$
The accepted step is the largest power $\beta^k \le \bar\gamma_c$. Because $\bar\gamma_c \ge 2(1-c)/C$, every accepted step satisfies
$$\beta \cdot \frac{2(1-c)}{C} < \gamma_n \le \bar\gamma_c. \tag{16}$$

### Effect of $\beta$ (with $c = 0.1$)

- **Small $\beta$** (0.1, 0.2): the grid $\{1, \beta, \beta^2, \dots\}$ is coarse, so the accepted step often lands **far below** $\bar\gamma_c$. For $\beta = 0.1$ the step is usually $10^{-3} = 0.565/C$, which is smaller than FSS with $p = 1$. Few reductions, but slow convergence (gap $0.135$).
- **$\beta = 0.5$**: a finer grid that lands within a factor of 2 of $\bar\gamma_c$. Steps of $1.1/C$ to $8.8/C$ were observed, which gives the **best gap ($1.13\times10^{-4}$)**.
- **$\beta = 0.9$**: very fine grid, so steps land close to $\bar\gamma_c$ (about $1.7/C$). Starting from $\gamma = 1$, though, it needs about 55 reductions per iteration, because $0.9^{60} \approx 1/C$. Many evaluations, and the gap is $0.063$.
- **Long steps:** BLS sometimes accepts $\gamma > 2/C$, beyond the fixed-step stability limit, for example $4.4/C$ or $8.8/C$. This is safe because Armijo checks decrease **along the current $g$** only. Occasional long steps make large progress along small-eigenvalue directions (similar to the Barzilai–Borwein effect), which explains why BLS with $\beta = 0.5$ beats ELS here.

### Effect of $c$ (with $\beta = 0.5$)

- A larger $c$ makes the condition stricter: (15) gives $\bar\gamma_c \propto (1 - c)$, so going from $c = 0.01$ to $c = 0.1$ shrinks the threshold by a factor of $0.9/0.99 \approx 0.91$.
- **When results can be identical:** if no power of $\beta$ falls between $\bar\gamma_{0.1}$ and $\bar\gamma_{0.01}$, both settings accept the same step. Here the steps are **not** identical. The trajectories diverge after the first iterations, with final gaps of $6.24\times10^{-4}$ for $c = 0.01$ and $1.13\times10^{-4}$ for $c = 0.1$. A stricter $c$ does not necessarily mean slower convergence; the difference comes from which grid points get selected.

### Cost (Table 1)

- The evaluation count follows an exact identity, which is worth stating in the report:
  $$\#\text{evals} = \underbrace{1}_{L(x_0)} + \underbrace{30}_{\text{accepted trials (cached as } L(x_n))} + \#\text{reductions}. \tag{17}$$
- Runtime scales roughly with the number of evaluations ($\beta = 0.9$: 1658 evaluations, about 10.9 ms; others: 116–269 evaluations, about 1–2 ms).
- Differences between $\beta = 0.1$ and $0.2$ (about 1.0–1.2 ms) are within the IQR, so do not claim a meaningful difference there.

---

## Problem 6: Comparison and recommendation

| Method | $L(x_{30}) - L_{\text{opt}}$ | Updates | Obj. evals | Median runtime |
|---|---|---|---|---|
| FSS, $p = 1$ | $9.81 \times 10^{-1}$ | 30 | 0 | $\approx 0.15$ ms |
| ELS | $4.44 \times 10^{-2}$ | 30 | 0 | $\approx 0.26$ ms |
| BLS, $c = 0.1,\ \beta = 0.5$ | $1.13 \times 10^{-4}$ | 30 | 267 | $\approx 1.9$ ms |

**Work not captured by evaluation counts.** Matrix-vector products per iteration, with $K = 50$ and $N = 100$:
$$\text{FSS: } 1\ (K\!\times\!K), \qquad \text{ELS: } 2\ (K\!\times\!K), \qquad \text{BLS: } 1\ (K\!\times\!K) + (1 + \#\text{backtracks}_n)\ (N\!\times\!K). \tag{18}$$
- ELS performs zero objective evaluations but still pays for one extra $G g$ product per step.
- All runtimes **exclude** CSV loading, the precomputation of $G$, $v$, $C$, and the SVD benchmark. State this explicitly in the report.

**Log-scale plotting.** Gaps are clipped to $\varepsilon \cdot L_{\text{opt}}$ for plotting only. No values were actually clipped, and the tables use the raw values.

**Recommendation (to argue).**
- **Progress per iteration:** BLS > ELS > FSS.
- **Progress per unit time:** BLS costs about 12× FSS in runtime but reaches a gap about 8700× smaller. ELS is cheap and gains about 22× over FSS.
- Reasonable choices are **BLS (0.1, 0.5)** if the target is accuracy within 30 iterations, or **ELS** if per-iteration cost matters most. For an LS problem, ELS needs no tuning parameters.

**Limitation (pick one).**
- ELS's closed form exists only for quadratic objectives.
- BLS's advantage here depends on lucky long steps from the $\beta$-grid and the $\gamma = 1$ reset, and it may not transfer to other data.
- With $\kappa \approx 22$ this problem is well-conditioned, and the rankings can change for ill-conditioned or rank-deficient $A$.
- Microsecond-level timings in R depend on the platform and BLAS library.

---

## Report-format checklist

- [ ] Number every referenced equation; caption every figure with its **main conclusion**.
- [ ] Fig. 1–6 have labeled axes, legends, and the $L_{\text{opt}}$ or $1/C$ dashed benchmarks.
- [ ] State the tolerance ($10^{-8}$), the rank tolerance, the timing method (1 warm-up + 21 timed runs, median and IQR), and the software versions.
- [ ] Contribution statement (if working in a team), references (including lecture slides), and a code appendix.
- [ ] Resolve the Python-vs-R requirement.

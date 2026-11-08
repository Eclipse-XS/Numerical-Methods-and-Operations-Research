# Cluster 02.06 — Nonlinear Unconstrained Optimization

## 1. Computational Problem
Minimizing smooth multivariable objective functions $\min_{\mathbf{x} \in \mathbb{R}^n} f(\mathbf{x})$ using zero-order heuristic search, first-order gradient search, and quasi-Newton second-order methods.

## 2. Algorithms & Methods
- **Nelder-Mead Regular Simplex (`Nelder_Mead/`):** Direct geometric search in $\mathbb{R}^n$ using reflection, expansion, and contraction of a simplex polytope without derivative evaluation.
- **Fletcher-Reeves Conjugate Gradient (`Fletcher_Reeves/`):** First-order gradient optimization generating mutually conjugate search directions $\mathbf{p}_k = -\mathbf{g}_k + \beta_k \mathbf{p}_{k-1}$ where $\beta_k = \frac{\|\mathbf{g}_k\|^2}{\|\mathbf{g}_{k-1}\|^2}$.
- **Davidon-Fletcher-Powell (DFP) (`DFP/`):** Quasi-Newton optimization updating an inverse Hessian approximation matrix $H_{k+1} = H_k + \Delta H_k$ satisfying the quasi-Newton secant condition.
- **Parabolic Interpolation (`Other_Methods/`):** One-dimensional line search finding the extremum of a quadratic polynomial passing through three bracketed points.

## 3. Included Projects & Subdirectories

| Subdirectory | Project | Algorithm | Distinctive Features | Assessment |
| :--- | :--- | :--- | :--- | :--- |
| `Nelder_Mead/` | `DO_RGR` | Nelder-Mead Regular Simplex | 2D geometric simplex search with automated reflection | **Canonical Reference (Zero-order)** |
| `Fletcher_Reeves/` | `DO_RGR4` | Fletcher-Reeves Gradient | Analytical gradient via SymPy, exact line minimization | **Canonical Reference (First-order)** |
| `DFP/` | `RGR4_Maks` | Davidon-Fletcher-Powell (DFP) | Quasi-Newton inverse Hessian update matrix trace | **Canonical Reference (Quasi-Newton)** |
| `Other_Methods/` | `rgr_21` | Parabolic Line Search | Parabolic minimum localization for Variant 21 | Historical Variant |
| `Other_Methods/` | `Operations_Research_RGR_Variant_21` | Line Search Substitution Trace | Step-by-step substitution trace for Variant 21 | Historical Deliverable |

## 4. Reference Implementations & Algorithmic Notes
Each subdirectory preserves a distinct mathematical optimization paradigm. `DO_RGR` provides derivative-free search; `DO_RGR4` leverages analytical gradient information; and `RGR4_Maks` models quadratic convergence through curvature estimation.

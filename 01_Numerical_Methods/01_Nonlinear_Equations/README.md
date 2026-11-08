# Cluster 01.01 — Nonlinear Equations (Root Finding)

## 1. Computational Problem
Finding real roots $x^*$ such that $f(x) = 0$ for single-variable nonlinear and transcendental algebraic equations over a bounded interval $[a, b]$ with specified convergence tolerance $\epsilon$.

## 2. Algorithms & Methods
- **Simple Fixed-Point Iteration:** Transforming $f(x) = 0$ into $x = \varphi(x)$ where $|\varphi'(x)| < 1$.
- **Newton-Raphson Method:** Tangent approximation $x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)}$ with quadratic convergence.
- **Secant / Chord Method:** Finite-difference secant approximation eliminating the requirement for analytical derivatives.
- **Symbolic Calculus:** Automated analytical differentiation using SymPy (`diff(f, x)`).

## 3. Included Projects

| Project | Algorithm | Distinctive Features | Relationship | Assessment |
| :--- | :--- | :--- | :--- | :--- |
| `Function CHM1` | Symbolic Newton & Chord | Object-oriented `Argument` class, analytical derivatives via SymPy, step table | Primary Reference | **Canonical Reference** (strongest OOP design) |
| `Non-linear-equations` | Newton, Secant, Iteration | Side-by-side benchmark of 3 methods with $\epsilon = 10^{-4}$ tolerance | Comparative Benchmark | High Academic Value |
| `CHM1` | Procedural Iteration | Basic interactive console input loop with manual evaluation | Earlier Draft | Historical / Redundant |

## 4. Reference Implementation & Unique Features
- **Strongest Reference:** `Function CHM1`. It decouples function evaluation, interval bounds, and step-by-step table generation using structured object modeling.
- **Unique Features:** `Non-linear-equations` provides a direct convergence comparison between three classical iterative techniques on $f(x) = 2\ln(x) - \frac{1}{x}$.

## 5. Known Limitations & Correctness
- `CHM1` lacks input validation and may enter an infinite loop if initial bracket $[a, b]$ violates the intermediate value theorem ($f(a) \cdot f(b) > 0$).

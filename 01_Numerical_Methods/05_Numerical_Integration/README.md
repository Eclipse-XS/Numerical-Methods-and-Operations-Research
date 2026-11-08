# Cluster 01.05 — Numerical Quadrature & Integration

## 1. Computational Problem
Approximating definite integrals $I = \int_a^b f(x) \, dx$ where analytical antiderivatives are complex or unavailable.

## 2. Algorithms & Methods
- **Composite Rectangle Rules:** Left, Right, and Midpoint approximations.
- **Composite Trapezoidal Rule:** Piecewise linear approximation.
- **Composite Simpson's 1/3 Rule:** Piecewise parabolic quadrature with $O(h^4)$ error term.

## 3. Included Projects

| Project | Algorithm | Distinctive Features | Relationship | Assessment |
| :--- | :--- | :--- | :--- | :--- |
| `6CHM` | Composite Quadrature Benchmark | Unified execution of 5 quadrature rules with error estimation | Standalone Lab | **Canonical Reference** |

## 4. Reference Implementation & Unique Features
- Provides direct numerical comparisons across varying grid step counts $N$, demonstrating the superior rate of convergence of Simpson's rule.

## 5. Known Limitations
- Grid step $N$ is fixed within the script; Runge's rule for adaptive step-size error control is not automated.

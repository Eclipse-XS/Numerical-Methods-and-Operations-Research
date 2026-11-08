# Cluster 01.04 — Systems of Nonlinear Equations

## 1. Computational Problem
Simultaneously solving a vector system of $m$ non-linear equations $\mathbf{F}(\mathbf{x}) = \mathbf{0}$ where $\mathbf{x} = (x_1, x_2)^T$.

## 2. Algorithms & Methods
- **Multidimensional Newton-Raphson:**
  $$\mathbf{x}_{k+1} = \mathbf{x}_k - J^{-1}(\mathbf{x}_k) \mathbf{F}(\mathbf{x}_k)$$
- **Jacobian Matrix Computation:** Exact analytical computation of partial derivatives $\frac{\partial F_i}{\partial x_j}$ using SymPy.

## 3. Included Projects

| Project | Algorithm | Distinctive Features | Relationship | Assessment |
| :--- | :--- | :--- | :--- | :--- |
| `5_CHM` | 2D Newton-Raphson | Analytical Jacobian matrix generation via SymPy, Tabulate iteration report | Standalone Lab | **Canonical Reference** |

## 4. Reference Implementation & Unique Features
- Dynamic calculation of the $2 \times 2$ Jacobian determinant and matrix inversion at each iteration step, outputting full numerical convergence telemetry.

## 5. Known Limitations
- Convergence is locally quadratic and depends strictly on initial point selection $(\mathbf{x}_0)$. Singular Jacobian states are not handled with damping.

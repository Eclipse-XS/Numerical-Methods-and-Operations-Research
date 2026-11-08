# Cluster 01.07 — Ordinary Differential Equations (ODEs)

## 1. Computational Problem
Solving initial value problems (IVP) $y' = f(x, y)$ and two-point boundary value problems (BVP) for second-order differential equations.

## 2. Algorithms & Methods
- **Runge-Kutta 4th Order (RK4):** Classical 4th-order multi-stage integrator with $O(h^4)$ local truncation error.
- **Finite-Difference Method (BVP):** Discretizing $y'' + p(x)y' + q(x)y = f(x)$ on a uniform mesh using central difference formulas, solving the resulting tridiagonal SLAE.

## 3. Included Projects

| Project | Algorithm | Distinctive Features | Relationship | Assessment |
| :--- | :--- | :--- | :--- | :--- |
| `4_RGR_CHM` | Runge-Kutta 4th Order | Step $h=0.01$, interval $[0, 2]$, Matplotlib trajectory plot | Term Project (IVP) | **Canonical Reference (IVP)** |
| `5_RGR_CHM` | Finite Differences (BVP) | Step $h=0.25$, interval $[0, 1]$, boundary conditions, mesh plot | Term Project (BVP) | **Canonical Reference (BVP)** |

## 4. Reference Implementation & Unique Features
- Together, these two projects cover the dual foundations of numerical differential equations: IVP integration via RK4 and boundary value discretization via finite differences.

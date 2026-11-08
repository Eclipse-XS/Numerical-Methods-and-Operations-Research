# Cluster 01.06 — Function Approximation & Spline Interpolation

## 1. Computational Problem
Constructing smooth piecewise interpolating polynomials passing through a discrete set of coordinate nodes $(x_i, y_i)_{i=0}^n$.

## 2. Algorithms & Methods
- **Cubic Spline Interpolation:** Piecewise $C^2$-continuous cubic polynomials $S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3$ matching 1st and 2nd derivatives at internal knots.
- **Boundary Conditions:** Natural splines ($S''(x_0) = S''(x_n) = 0$).

## 3. Included Projects

| Project | Algorithm | Distinctive Features | Relationship | Assessment |
| :--- | :--- | :--- | :--- | :--- |
| `7CHM` | Cubic Splines (SciPy) | High-resolution Matplotlib curve comparison, node tabularization | Primary Reference | **Canonical Reference** |
| `rozum_7chm` | Cubic Splines & SymPy | Symbolic polynomial integration and LaTeX equation rendering | Peer Variant | High Analytical Value |

## 4. Reference Implementation & Unique Features
- `7CHM` provides professional visual verification using SciPy and Matplotlib.
- `rozum_7chm` contains unique analytical integration steps using SymPy to compute the exact area under the spline.

## 5. Known Limitations
- Both implementations require strictly sorted knot arrays along $x$.

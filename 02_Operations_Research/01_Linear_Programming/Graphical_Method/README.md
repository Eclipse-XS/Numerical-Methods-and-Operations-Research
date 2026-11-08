# Cluster 02.01.A — Linear Programming: 2D Graphical Method

## 1. Computational Problem
Solving two-variable linear programs $\max / \min \{ c_1 x_1 + c_2 x_2 \}$ subject to linear half-plane inequalities $A \mathbf{x} \le \mathbf{b}, \mathbf{x} \ge 0$.

## 2. Algorithms & Methods
- **Half-Plane Intersection:** Pairwise line intersections via `np.linalg.solve`.
- **Feasible Polygon Filtering:** Checking candidate vertices against all inequality constraints.
- **Gradient Level Lines:** Evaluating $Z(\mathbf{x})$ on vertices and tracing objective normal vector $\mathbf{c}$.

## 3. Included Projects

| Project | Algorithm | Distinctive Features | Relationship | Assessment |
| :--- | :--- | :--- | :--- | :--- |
| `DO_1` | 2D Graphical Solver | Author signature (OI-23), clean polygon validation, Matplotlib plot | Primary Reference | **Canonical Reference** |
| `DO1_Ostap` | 2D Graphical Solver | Formatted terminal output tailored for Ostap | Parameter Variant | Historical Variant |
| `DO_Denys1` | 2D Graphical Solver | Variant equations ($2x - 3y = 6$) for Denys | Parameter Variant | Historical Variant |
| `DO_Denys1_GORLOVOY` | 2D Graphical Solver | Grid-sampling verification mesh for Denys | Parameter Variant | Historical Variant |

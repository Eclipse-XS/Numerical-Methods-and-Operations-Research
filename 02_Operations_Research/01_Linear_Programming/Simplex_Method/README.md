# Cluster 02.01.B — Linear Programming: Simplex Method

## 1. Computational Problem
Solving canonical general linear programs $\max \mathbf{c}^T \mathbf{x}$ subject to $A \mathbf{x} = \mathbf{b}, \mathbf{x} \ge \mathbf{0}$ with non-negative right-hand sides.

## 2. Algorithms & Methods
- **Tableau Transformation:** Maintaining basic and non-basic variable partitions.
- **Optimality Criterion:** Checking reduced cost row $\Delta_j = z_j - c_j \ge 0$.
- **Pivot Selection:** Dantzig's rule for entering variable ($\min \Delta_j$) and minimum ratio test $\min_i \{ b_i / a_{ik} \mid a_{ik} > 0 \}$ for leaving variable.

## 3. Included Projects

| Project | Algorithm | Distinctive Features | Relationship | Assessment |
| :--- | :--- | :--- | :--- | :--- |
| `DO_2` | Standard Simplex | Canonical reference tableau solver by author | Primary Reference | **Canonical Reference (Python)** |
| `do` | C++ Simplex Solver | Standalone C++ console implementation (`FileName.cpp`, 326 lines) | Unique Language | **Canonical Reference (C++)** |
| `DO_Another` | Simplex + Graphical | Hybrid solver featuring both algebraic simplex and 2D feasible plot | Extended Tool | High Educational Value |
| `Another_DO2_94` | Standard Simplex | Clean ASCII table renderer for Variant 94 | Parameter Variant | Historical Variant |
| `DO2_Ostap` | Standard Simplex | ASCII border formatting tailored for Ostap | Parameter Variant | Historical Variant |
| `DO2_Denys` | Standard Simplex | Colored terminal tableau output for Denys | Parameter Variant | Historical Variant |

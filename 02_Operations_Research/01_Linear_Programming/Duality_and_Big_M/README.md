# Cluster 02.01.C — Linear Programming: Duality & Big-M Method

## 1. Computational Problem
Solving linear programs requiring artificial starting bases and constructing mathematically equivalent dual problems.

## 2. Algorithms & Methods
- **Primal-Dual Transformation:** Transposing coefficient matrix $A^T$, inverting objective sense (max $\leftrightarrow$ min), and mapping constraints to dual variables.
- **Big-M Method:** Introducing artificial variables with penalty $-M$ (or $+M$) into the objective to force an initial identity basis.

## 3. Included Projects

| Project | Algorithm | Distinctive Features | Relationship | Assessment |
| :--- | :--- | :--- | :--- | :--- |
| `DO_3(final)` | Big-M + Duality | Symbolic $M$ penalty, complete canonical solver | Primary Reference | **Canonical Reference (Python)** |
| `DO_BIG-M` | C++ Big-M Solver | Standalone C++ console implementation of Big-M | Unique Language | **Canonical Reference (C++)** |
| `DO_3_Duality` | Duality Builder | Dedicated primal-to-dual matrix converter | Utility Tool | Modular Component |
| `DO_3` | Big-M Solver | Early procedural draft | Precursor | Historical Draft |
| `Big_M_DO3` | Big-M Solver | Alternative compact implementation | Variant | Historical Variant |
| `DO3_OStap` | Big-M Solver | Customized variant for Ostap | Parameter Variant | Historical Variant |
| `DO3_Denys` | Big-M Solver | Customized variant for Denys | Parameter Variant | Historical Variant |
| `DO_3_Denys` | Big-M Solver | Alternative variant for Denys | Parameter Variant | Historical Variant |
| `Another_DO35_3` | Big-M Solver | Customized variant for Variant 35 | Parameter Variant | Historical Variant |
| `Another_DO3_72` | Big-M Solver | Customized variant for Variant 72 | Parameter Variant | Historical Variant |
| `Another_DO3_94` | Big-M Solver | Customized variant for Variant 94 | Parameter Variant | Historical Variant |

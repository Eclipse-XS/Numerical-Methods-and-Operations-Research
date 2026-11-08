# Cluster 02.03 — Integer Linear Programming (Gomory's Cutting Planes)

## 1. Computational Problem
Solving linear programs with integer constraints $\mathbf{x} \in \mathbb{Z}_+^n$ using cutting-plane methods.

## 2. Algorithms & Methods
- **Continuous Relaxation:** Solving LP via dual simplex.
- **Gomory Fractional Cuts:** Generating valid cutting planes from rows with fractional basic solutions:
  $$\sum_{j \in N} \{ a_{ij} \} x_j \ge \{ b_i \}$$
- **Exact Symbolic Arithmetic:** Using SymPy exact rational representations (`Rational`) to eliminate floating-point truncation errors during cut derivation.

## 3. Included Projects

| Project | Algorithm | Distinctive Features | Relationship | Assessment |
| :--- | :--- | :--- | :--- | :--- |
| `THE_BEST_CALCULATOR_DO5_GOMORY_CUTTING` | Symbolic Gomory | Exact rational arithmetic, Word (`.docx`) table export | Master Solver | **Canonical Master Tool** |
| `DO_5` | Gomory Cutting Plane | Original student reference lab by author | Primary Reference | **Canonical Reference** |
| `NEW_GOMORY` | Symbolic Gomory | 98.9% near-duplicate of master calculator without docx | Variant | Redundant Tool |
| `CALCULATOR_DO5` | Symbolic Gomory | Interactive Gomory matrix solver | Prototype | Historical Tool |
| `DO5_Maks(3)` | Symbolic Gomory | Cleaned final deliverable for student Maks | Clean Deliverable | Historical Variant |
| `DO5_Maks` | Symbolic Gomory | Version with debug prints | Intermediate | Superseded by (3) |
| `DO5_Maks(2)` | Fractions Draft | Incomplete exploratory draft with `math.floor` | Broken Draft | **DELETE_CANDIDATE** |
| `DO5_Maks(4)` | Truncated Draft | Truncated script fragment | Broken Draft | **DELETE_CANDIDATE** |
| `DO5_Ostap` | Symbolic Gomory | Customized variant for Ostap | Parameter Variant | Historical Variant |
| `DO5_Denys` | Float Big-M Gomory | Early attempt using $M=10^{12}$ | Flawed Approach | Historical Variant |
| `DO5_Denys_Improved` | Symbolic Gomory | Improved symbolic rewrite for Denys | Refinement | Historical Variant |
| `DO5_Kremitets` | Symbolic Gomory | Customized variant for Kremitets | Parameter Variant | Historical Variant |
| `Another_72_DO5` | Symbolic Gomory | Customized variant for Variant 72 | Parameter Variant | Historical Variant |
| `Another_94_DO5` | Symbolic Gomory | Customized variant for Variant 94 | Parameter Variant | Historical Variant |
| `Another_DO5_19` | Symbolic Gomory | Customized variant for Variant 19 | Parameter Variant | Historical Variant |
| `Another_DO5_35` | Symbolic Gomory | Customized variant for Variant 35 | Parameter Variant | Historical Variant |
| `Another_DO5_94` | Symbolic Gomory | Customized variant for Variant 94 | Parameter Variant | Historical Variant |

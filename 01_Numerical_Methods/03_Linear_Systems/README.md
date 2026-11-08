# Cluster 01.03 — Systems of Linear Algebraic Equations (SLAE)

## 1. Computational Problem
Solving $A x = b$ for dense square coefficient matrices $A \in \mathbb{R}^{n \times n}$ and right-hand side vectors $b \in \mathbb{R}^n$.

## 2. Algorithms & Methods
- **Gauss-Jordan Elimination:** Systematic forward elimination and backward substitution transforming augmented matrix $[A \mid b]$ directly into reduced row echelon form $[I \mid x^*]$.
- **Pivoting:** Partial row pivoting to eliminate zero-division risks on diagonal entries.

## 3. Included Projects

| Project | Algorithm | Distinctive Features | Relationship | Assessment |
| :--- | :--- | :--- | :--- | :--- |
| `Gauss-Jordan` | Gauss-Jordan Elimination | Step-by-step augmented matrix pretty-printing, ANSI colored pivots | Standalone Lab | **Canonical Reference** |

## 4. Reference Implementation & Unique Features
- Prints every elementary row operation ($R_i \leftarrow R_i - c \cdot R_j$) with clear terminal formatting and verification.

## 5. Known Limitations
- ANSI escape codes are hardcoded without `colorama.init()`, which may display raw escape characters on legacy Windows terminals.

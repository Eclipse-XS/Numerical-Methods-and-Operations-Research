# Cluster 02.04.A — Fractional-Linear Programming (Algebraic Approach)

## 1. Computational Problem
Maximizing or minimizing a ratio of linear functions:
$$F(\mathbf{x}) = \frac{\mathbf{c}^T \mathbf{x} + \alpha}{\mathbf{d}^T \mathbf{x} + \beta} \quad \text{subject to } A \mathbf{x} \le \mathbf{b}, \mathbf{x} \ge \mathbf{0}$$

## 2. Algorithms & Methods
- **Charnes-Cooper Transformation:** Introducing transformation variables:
  $$t = \frac{1}{\mathbf{d}^T \mathbf{x} + \beta}, \quad \mathbf{y} = t \cdot \mathbf{x}$$
  Yielding an equivalent linear program in $(\mathbf{y}, t)$ space:
  $$\max \mathbf{c}^T \mathbf{y} + \alpha t \quad \text{s.t. } A \mathbf{y} - \mathbf{b} t \le \mathbf{0}, \; \mathbf{d}^T \mathbf{y} + \beta t = 1, \; \mathbf{y} \ge \mathbf{0}, t \ge 0$$

## 3. Included Projects

| Project | Algorithm | Distinctive Features | Relationship | Assessment |
| :--- | :--- | :--- | :--- | :--- |
| `DO6_simplex` | Charnes-Cooper Simplex | Symbolic transformation and complete canonical tableau solver | Primary Reference | **Canonical Reference** |
| `DO_6` | Fractional LP Solver | Original student reference lab by author | Primary Reference | High Value Reference |
| `DO6_Nastya` | Charnes-Cooper Simplex | Formatted output with symbolic expressions | Peer Variant | Clean Variant |
| `DO6_Maks` | Charnes-Cooper Simplex | Customized variant for Maks | Parameter Variant | Historical Variant |
| `DO6_Ostap` | Charnes-Cooper Simplex | Customized variant for Ostap | Parameter Variant | Historical Variant |
| `DO6_Denys` | Charnes-Cooper Simplex | Colored tableau variant for Denys | Parameter Variant | Historical Variant |
| `BUZZKILLER_DO6_35` | Charnes-Cooper Simplex | Customized variant for Variant 35 | Parameter Variant | Historical Variant |
| `Another_72_DO6` | Charnes-Cooper Simplex | Customized variant for Variant 72 | Parameter Variant | Historical Variant |
| `Another_94_DO6` | Charnes-Cooper Simplex | Customized variant for Variant 94 | Parameter Variant | Historical Variant |

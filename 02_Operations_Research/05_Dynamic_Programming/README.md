# Cluster 02.05 — Dynamic Programming (Bellman Allocation)

## 1. Computational Problem
Multi-stage optimal resource allocation: distributing capital $S$ across $N$ industrial enterprises to maximize total profit $\sum f_i(x_i)$ where $\sum x_i \le S$.

## 2. Algorithms & Methods
- **Bellman's Principle of Optimality:** Backward recursion:
  $$F_k(s) = \max_{0 \le x \le s} \{ f_k(x) + F_{k-1}(s - x) \}$$
- **State Tabularization:** Maintaining state-decision tables for each allocation stage.

## 3. Included Projects

| Project | Algorithm | Distinctive Features | Relationship | Assessment |
| :--- | :--- | :--- | :--- | :--- |
| `DO7_CALCIULATOR_FOR_Students` | Multi-Stage Bellman Solver | Generalized solver with Tabulate multi-stage tables | Master Solver | **Canonical Master Tool** |
| `DO_7` | Bellman Allocation | Original student reference lab by author | Primary Reference | **Canonical Reference** |
| `DO6_Nasty` | Multi-Stage Bellman Solver | Grounded analysis proves this is Lab 7 despite "DO6" name | Misnamed Variant | Historical Variant |
| `DO7_Ostap` | Bellman Allocation | Customized variant for Ostap | Parameter Variant | Historical Variant |
| `DO7_Denys` | Bellman Allocation | Customized variant for Denys | Parameter Variant | Historical Variant |

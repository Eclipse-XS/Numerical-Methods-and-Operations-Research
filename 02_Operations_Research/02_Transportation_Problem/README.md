# Cluster 02.02 — Transportation Problem & Method of Potentials

## 1. Computational Problem
Optimizing shipping distribution from $m$ suppliers with capacities $a_i$ to $n$ consumers with demands $b_j$ minimizing total cost $\sum c_{ij} x_{ij}$.

## 2. Algorithms & Methods
- **Balance Verification:** Balancing open models ($\sum a_i \ne \sum b_j$) by appending fictitious dummy suppliers or consumers.
- **Initial Basic Feasible Solution:** North-West Corner method and Minimum Cost element rule.
- **Optimality & Potentials:** Evaluating dual potentials $u_i + v_j = c_{ij}$ for basic cells and checking $\Delta_{ij} = c_{ij} - (u_i + v_j) \ge 0$.
- **Stepping-Stone Cycle:** Detecting closed loops to reallocate shipment volume $\theta$.

## 3. Included Projects

| Project | Algorithm | Distinctive Features | Relationship | Assessment |
| :--- | :--- | :--- | :--- | :--- |
| `DO_4_Maks` | Potentials & Balancing | Git-tracked repo, colored terminal display, balance checks | Primary Reference | **Canonical Reference** |
| `DO_4` | Potentials Method | Author's original reference lab | Precursor | High Value Reference |
| `Another_DO5_72` | Potentials Method | Grounded analysis proves this is Lab 4 despite "DO5" name | Misnamed Variant | Historical Variant |
| `Another_DO19_4` | Potentials Method | Customized variant for Variant 19 | Parameter Variant | Historical Variant |
| `DO4_Ostap` | Potentials Method | Customized variant for Ostap | Parameter Variant | Historical Variant |
| `DO4_Denys` | Potentials Method | Customized variant for Denys | Parameter Variant | Historical Variant |
| `Another_DO4_35` | Potentials Method | Customized variant for Variant 35 | Parameter Variant | Historical Variant |
| `Another_DO4_54` | Potentials Method | Customized variant for Variant 54 | Parameter Variant | Historical Variant |
| `Another_DO4_94` | Potentials Method | Customized variant for Variant 94 | Parameter Variant | Historical Variant |

# Cluster 01.02 — Polynomial Root Isolation (Sturm Chains)

## 1. Computational Problem
Determining the exact number of distinct real roots of a polynomial $P(x)$ in a given interval $(a, b]$ without computing the roots directly, avoiding root-clustering errors.

## 2. Algorithms & Methods
- **Sturm's Theorem:** Constructing a canonical Sturm chain:
  $$P_0(x) = P(x), \quad P_1(x) = P'(x), \quad P_{k}(x) = -\text{rem}(P_{k-2}, P_{k-1})$$
- **Sign Variations:** Counting sign changes $V(x)$ at interval boundaries to calculate exact real roots $N(a, b) = V(a) - V(b)$.

## 3. Included Projects

| Project | Algorithm | Distinctive Features | Relationship | Assessment |
| :--- | :--- | :--- | :--- | :--- |
| `Sturm Method` | Sturm Chain Construction | SymPy polynomial division `Poly.div`, automated chain evaluation, author signature (OI-23) | Standalone Lab | **Canonical Reference** |

## 4. Reference Implementation & Unique Features
- Clean, mathematically exact implementation using SymPy's algebraic polynomial engine. Formats intermediate polynomials and sign variation vectors explicitly.

## 5. Known Limitations
- Hardcoded polynomial definition inside source code. Requires manual script edits to test other polynomials.

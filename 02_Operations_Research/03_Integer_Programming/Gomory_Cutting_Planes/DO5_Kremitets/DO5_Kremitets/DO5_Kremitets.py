import sympy as sp
import numpy as np
M = sp.Symbol("M")
def setup_lp(c, A, b, signs, goal='max'):
    n, A, b = len(c), sp.Matrix(A), sp.Matrix(b).reshape(len(b), 1)
    names, basis, art_vars, tab = [f"x{i+1}" for i in range(n)], [], [], A
    idx = n + 1

    for i, s in enumerate(signs):
        if s == "<=":
            tab = tab.row_join(sp.Matrix.eye(len(b))[:, i])
            names.append(f"x{idx}"); basis.append(f"x{idx}"); idx += 1
        elif s == ">=":
            tab = tab.row_join(-sp.Matrix.eye(len(b))[:, i])
            names.append(f"x{idx}"); idx += 1
            tab = tab.row_join(sp.Matrix.eye(len(b))[:, i])
            names.append(f"w{len([v for v in names if v.startswith('w')]) + 1}"); basis.append(names[-1])
            art_vars.append((i, tab.cols - 1))
        elif s == "=":
            tab = tab.row_join(sp.Matrix.eye(len(b))[:, i])
            names.append(f"w{len([v for v in names if v.startswith('w')]) + 1}"); basis.append(names[-1])
            art_vars.append((i, tab.cols - 1))

    tab = tab.row_join(b)
    cost = sp.Matrix([0]*(tab.cols - 1) + [0])
    for i in range(n): cost[i] = sp.Rational(c[i])

    m = 1 if goal == 'min' else -1
    for i, _ in art_vars:
        row = tab.row(i)
        for j in range(len(cost)-1):
            cost[j] += m * M * (-row[j])
        cost[-1] += m * M * row[-1]

    for _, col in art_vars: cost[col] = 0
    cost = cost.applyfunc(lambda x: -x)  

    return tab, cost, names, basis
def simplex(tab, c, names, basis, goal='max'):
    tab = np.array(tab.tolist(), dtype=object)
    c = [sp.sympify(x) for x in c]
    m, n = tab.shape[0], tab.shape[1] - 1

    def pivot_col(c, goal='max'):
        vals = [float(sp.sympify(x).subs(M, 1000)) for x in c[:-1]]
        if goal == 'max':
            return vals.index(min(vals)) if min(vals) < 0 else None
        return max((i for i, v in enumerate(vals) if v > 0), key=lambda i: vals[i], default=None)

    while (col := pivot_col(c, goal)) is not None:
        ratios = [(tab[i, -1]/tab[i, col] if tab[i, col] > 0 else float('inf')) for i in range(m)]
        if all(r == float('inf') for r in ratios): return None, None, None, None

        row = min((i for i, r in enumerate(ratios)), key=lambda i: ratios[i])
        pivot = tab[row, col]
        tab[row] = [sp.simplify(x/pivot) for x in tab[row]]
        
        for i in range(m):
            if i != row:
                factor = tab[i, col]
                tab[i] = [sp.simplify(tab[i][j] - factor * tab[row][j]) for j in range(tab.shape[1])]
        
        factor = c[col]
        c = [sp.simplify(c[j] - factor * tab[row][j]) for j in range(len(c))]
        basis[row] = names[col]

    sol = [0]*n
    for i in range(m):
        if basis[i] in names[:n]:
            sol[names.index(basis[i])] = tab[i, -1]
    return tab, c, basis, sol
def dual_simplex(tab, c, names, basis):
    m, n = tab.shape[0], tab.shape[1] - 1

    while any(tab[i, -1] < 0 for i in range(m)):
        row = min((i for i in range(m) if tab[i, -1] < 0), key=lambda i: tab[i, -1])
        candidates = [(abs(c[j]/tab[row, j]), j) for j in range(n) if tab[row, j] < 0]

        if not candidates: return None, None, None, None
        _, col = min(candidates, key=lambda x: float(x[0].subs(M, 1000) if isinstance(x[0], sp.Basic) else x[0]))
        
        pivot = tab[row, col]
        tab[row] = [sp.simplify(x/pivot) for x in tab[row]]
        for i in range(m):
            if i != row:
                factor = tab[i, col]
                tab[i] = [sp.simplify(tab[i][j] - factor * tab[row][j]) for j in range(tab.shape[1])]
        
        factor = c[col]
        c = [sp.simplify(c[j] - factor * tab[row][j]) for j in range(len(c))]
        basis[row] = names[col]

    sol = [0]*n
    for i in range(m):
        if basis[i] in names[:n]:
            sol[names.index(basis[i])] = tab[i, -1]
    return tab, c, basis, sol
def int_check(x): return np.isclose(float(x), round(float(x)))
def frac(x): return float(x) - np.floor(float(x))
def cutting(tab, c, names, basis):
    def find_cut_row(tab, names, basis):
        rows = [(i, frac(tab[i, -1])) for i, var in enumerate(basis) if var in ("x1", "x2") and not int_check(tab[i, -1])]
        return max(rows, key=lambda r: r[1], default=(None, None))[0]
    r = find_cut_row(tab, names, basis)
    if r is None: return tab, c, names, basis, None

    row = tab[r]
    cut = [-sp.nsimplify(frac(row[i])) for i in range(len(row)-1)]
    rhs = -sp.nsimplify(frac(row[-1]))

    new_var = f"x{len([v for v in names if v.startswith('x')]) + 1}"
    names.append(new_var); basis.append(new_var)

    new_row = cut + [1, rhs]
    tab = np.hstack((tab[:, :-1], np.zeros((tab.shape[0], 1)), tab[:, -1:]))
    tab = np.vstack((tab, np.array(new_row, dtype=object)))
    c = np.append(c[:-1], [0, c[-1]])

    terms = " + ".join([f"{sp.nsimplify(a)}*{b}" for a, b in zip(cut, names[:-1]) if a != 0])
    terms += f" + 1*{new_var} = {rhs}"

    return tab, c, names, basis, terms

A = [[1,2], [2,4], [6,1], [1,3]]
b = [5, 16, 6, 9]
c = [3,1]
signs = ['>=', '<=', '>=', '<=']
goal = 'min'
tab, z, names, basis = setup_lp(c, A, b, signs, goal)
tab, z, basis, sol = simplex(tab, z, names, basis, goal)
print("\nSimplex result:")
print(", ".join(f"{names[i]} = {sp.nsimplify(sol[i], rational=True)}" for i in range(len(c))))
F = sum(sp.Rational(c[i])*sol[i] for i in range(len(c)))
print(f"F(x) = {sp.nsimplify(F, rational=True)}")
if all(int_check(v) for v in sol):
    print("\nConclusion: All integer.")
else:
    print("\nIntermediate solution not integer, adding Gomory cuts...")
    while True:
        tab, z, names, basis, cut = cutting(tab, z, names, basis)
        if cut: print("\nAdded constraint:\n", cut)
        tab, z, basis, sol = dual_simplex(tab, z, names, basis)
        print("\nDual simplex result:")
        print(", ".join(f"{names[i]} = {sp.nsimplify(sol[i], rational=True)}" for i in range(len(c))))
        F = sum(sp.Rational(c[i])*sol[i] for i in range(len(c)))
        print(f"F(x) = {sp.nsimplify(F, rational=True)}")
        if all(int_check(v) for v in sol):
            print("\nConclusion: All integer.")
            break
        else:
            print("\nStill not integer, adding another cut...")
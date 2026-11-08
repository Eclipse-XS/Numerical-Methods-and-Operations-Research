import sympy as sp
import numpy as np
from sympy import Matrix, nsimplify, simplify, Rational, Symbol
from fractions import Fraction
from tabulate import tabulate

M = sp.Symbol("M")

def is_int(val):
    return np.isclose(float(val), round(float(val)))

def get_frac(val):
    return float(val) - np.floor(float(val))

def simplex_method(table, cost, var_names, basis, goal='max'):
    table = np.array(table.tolist(), dtype=object)
    cost = [sp.sympify(x) for x in cost]
    rows, cols = table.shape[0], table.shape[1] - 1

    while True:
        z_vals = [float(sp.sympify(expr).subs(M, 1000)) for expr in cost[:-1]]
        if goal == 'max':
            if min(z_vals) >= 0:
                break
            pivot_col = z_vals.index(min(z_vals))
        else:
            pos = [i for i, val in enumerate(z_vals) if val > 0]
            if not pos:
                break
            pivot_col = max(pos, key=lambda i: z_vals[i])

        ratios = []
        for i in range(rows):
            coeff, rhs = table[i, pivot_col], table[i, -1]
            ratio = rhs / coeff if coeff > 0 else float('inf')
            ratios.append(ratio)

        if all(r == float('inf') for r in ratios):
            return None, None, None, None

        pivot_row = ratios.index(min(ratios, key=float))
        pivot_elem = table[pivot_row, pivot_col]
        table[pivot_row] = [simplify(x / pivot_elem) for x in table[pivot_row]]

        for i in range(rows):
            if i != pivot_row:
                factor = table[i, pivot_col]
                table[i] = [simplify(table[i][j] - factor * table[pivot_row][j]) for j in range(table.shape[1])]

        factor = cost[pivot_col]
        cost = [simplify(cost[j] - factor * table[pivot_row][j]) for j in range(len(cost))]
        basis[pivot_row] = var_names[pivot_col]

    solution = [0] * cols
    for i in range(rows):
        if basis[i] in var_names[:cols]:
            idx = var_names.index(basis[i])
            solution[idx] = table[i, -1]

    return table, cost, basis, solution

def add_gomory_cut(table, cost, var_names, basis):
    candidates = []
    for i, var in enumerate(basis):
        if var in ["x1", "x2"]:
            val = table[i, -1]
            if not is_int(val):
                candidates.append((i, get_frac(val)))

    if not candidates:
        return table, cost, var_names, basis, False

    row_idx, _ = max(candidates, key=lambda x: x[1])
    row = table[row_idx]
    cut_coeff = [Fraction(-get_frac(row[i])).limit_denominator(1000) for i in range(len(row) - 1)]
    cut_row = cut_coeff + [1, Fraction(-get_frac(row[-1])).limit_denominator(1000)]

    table = np.hstack((table[:, :-1], np.zeros((table.shape[0], 1)), table[:, -1:]))
    table = np.vstack((table, np.array(cut_row, dtype=object)))
    cost = np.append(cost[:-1], [0, cost[-1]])

    cut_name = f"cut{sum(1 for v in var_names if v.startswith('cut')) + 1}"
    var_names.append(cut_name)
    basis.append(cut_name)

    return table, cost, var_names, basis, True

def dual_simplex_method(table, cost, var_names, basis):
    rows, cols = table.shape[0], table.shape[1] - 1
    while True:
        if all(table[i, -1] >= 0 for i in range(rows)):
            break
        pivot_row = min((i for i in range(rows) if table[i, -1] < 0), key=lambda i: table[i, -1])
        row = table[pivot_row, :-1]
        options = [(abs(cost[j] / row[j]), j) for j in range(cols) if row[j] < 0]
        if not options:
            return None, None, None, None
        _, pivot_col = min(options, key=lambda x: float(x[0].subs(M, 1000) if isinstance(x[0], sp.Basic) else x[0]))
        pivot_elem = table[pivot_row, pivot_col]
        table[pivot_row] = [simplify(x / pivot_elem) for x in table[pivot_row]]

        for i in range(rows):
            if i != pivot_row:
                factor = table[i, pivot_col]
                table[i] = [simplify(table[i][j] - factor * table[pivot_row][j]) for j in range(table.shape[1])]

        factor = cost[pivot_col]
        cost = [simplify(cost[j] - factor * table[pivot_row][j]) for j in range(len(cost))]
        basis[pivot_row] = var_names[pivot_col]

    solution = [0] * cols
    for i in range(rows):
        if basis[i] in var_names[:cols]:
            idx = var_names.index(basis[i])
            solution[idx] = table[i, -1]

    return table, cost, basis, solution

def build_table(c, A, b, signs, goal='max'):
    n = len(c)
    A, b = Matrix(A), Matrix(b).reshape(len(b), 1)
    table = A
    var_names = [f"x{i+1}" for i in range(n)]
    basis = []
    artificial = []
    next_idx = n + 1

    for i, sign in enumerate(signs):
        if sign in ("<=", ">="):
            slack = Matrix.zeros(len(b), 1)
            slack[i] = 1 if sign == "<=" else -1
            table = table.row_join(slack)
            var_names.append(f"x{next_idx}")
            if sign == "<=" : basis.append(f"x{next_idx}")
            next_idx += 1
        if sign in (">=", "="):
            art = Matrix.zeros(len(b), 1)
            art[i] = 1
            table = table.row_join(art)
            name = f"w{len([v for v in var_names if v.startswith('w')]) + 1}"
            var_names.append(name)
            basis.append(name)
            artificial.append((i, table.cols - 1))

    table = table.row_join(b)

    cost_row = Matrix([0] * (table.cols - 1) + [0])
    for i in range(len(c)):
        cost_row[i] = Rational(c[i])

    mult = 1 if goal == 'min' else -1
    for i, _ in artificial:
        row = table.row(i)
        for j in range(len(cost_row) - 1):
            cost_row[j] += mult * M * (-row[j])
        cost_row[-1] += mult * M * row[-1]

    for _, col in artificial:
        cost_row[col] = 0

    for i in range(len(cost_row) - 1):
        cost_row[i] *= -1

    return table, cost_row, var_names, basis

A = [[-1, 1], [6, 7], [2, -3]]
b = [3, 42, 6]
c = [4, -2]
signs = ['<=', '<=', '<=']
goal = 'max'

T, Z, names, base = build_table(c, A, b, signs, goal)
T, Z, base, sol = simplex_method(T, Z, names, base, goal)

print("\nРезультат звичайного симплекс-методу:")
sol_str = [(names[i], nsimplify(sol[i], rational=True)) for i in range(len(names)) if names[i] in ["x1", "x2"]]
print(tabulate(sol_str, headers=["Змінна", "Значення"], tablefmt="grid"))
F = nsimplify(sum(Rational(c[i]) * sol[i] for i in range(len(c))), rational=True)
print(f"F(x) = {F}\n")

cuts = 0
if not all(is_int(sol[i]) for i in range(len(c))):
    print("Метод Гоморі: починаємо додавати зрізи.")

while not all(is_int(sol[i]) for i in range(len(c))):
    cuts += 1
    print(f"\nЗріз {cuts}:")
    T, Z, names, base, added = add_gomory_cut(T, Z, names, base)
    if not added:
        print(" Не вдалося додати новий зріз. Завершення.")
        break
    T, Z, base, sol = dual_simplex_method(T, Z, names, base)
    sol_str = [(names[i], nsimplify(sol[i], rational=True)) for i in range(len(names)) if names[i] in ["x1", "x2"]]
    print(tabulate(sol_str, headers=["Змінна", "Значення"], tablefmt="grid"))
    F = nsimplify(sum(Rational(c[i]) * sol[i] for i in range(len(c))), rational=True)
    print(f"F(x) = {F}")

if cuts == 0:
    print("Метод Гоморі не застосовувався — розв’язок був цілочисельним.")
else:
    print("\nФінальний результат після методу Гоморі:")
    sol_str = [(names[i], nsimplify(sol[i], rational=True)) for i in range(len(names)) if names[i] in ["x1", "x2"]]
    print(tabulate(sol_str, headers=["Змінна", "Значення"], tablefmt="grid"))
    print(f"F(x) = {F}")
    print(f"Кількість зрізів Гоморі: {cuts}")

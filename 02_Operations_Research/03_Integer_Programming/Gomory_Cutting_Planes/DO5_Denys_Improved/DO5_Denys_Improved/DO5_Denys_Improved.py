import sympy as sp
import numpy as np
from sympy import Matrix, simplify, nsimplify
from tabulate import tabulate

M = sp.Symbol("M")

def build_canonical_table(c, A, b, signs, goal='max'):
    num_vars = len(c)
    A = Matrix(A)
    b = Matrix(b).reshape(len(b), 1)
    table = A.row_join(Matrix.eye(len(b)) * 0)
    var_names = [f"x{i+1}" for i in range(num_vars)]
    slack_count = 0
    artificial_count = 0
    artificial_vars = []
    basis = []

    for i, sign in enumerate(signs):
        if sign == "<=":
            slack_var = f"x{len(var_names)+1}"
            table[i, num_vars + slack_count] = 1
            var_names.append(slack_var)
            basis.append(slack_var)
            slack_count += 1
        elif sign == ">=":
            surplus_var = f"x{len(var_names)+1}"
            art_var = f"w{artificial_count + 1}"
            table[i, num_vars + slack_count] = -1
            table = table.row_join(Matrix.zeros(table.rows, 1))
            table[i, table.cols - 1] = 1
            var_names += [surplus_var, art_var]
            basis.append(art_var)
            slack_count += 1
            artificial_count += 1
            artificial_vars.append((i, art_var))
        elif sign == "=":
            art_var = f"w{artificial_count + 1}"
            table = table.row_join(Matrix.zeros(table.rows, 1))
            table[i, table.cols - 1] = 1
            var_names.append(art_var)
            basis.append(art_var)
            artificial_count += 1
            artificial_vars.append((i, art_var))

    table = table.row_join(b)
    cost_row = [0] * (table.cols - 1)
    for i in range(len(c)):
        cost_row[i] = sp.Rational(c[i])
    cost_row = Matrix(cost_row + [0])
    sign = 1 if goal == 'min' else -1
    for idx, a_var in artificial_vars:
        rhs = table[idx, -1]
        for j in range(table.cols - 1):
            if var_names[j] != a_var:
                cost_row[j] += sign * M * (-table[idx, j])
        cost_row[var_names.index(a_var)] = 0
        cost_row[-1] += sign * M * rhs
    for i in range(len(cost_row) - 1):
        cost_row[i] *= -1

    return table, cost_row, var_names, basis

def print_table_with_tabulate(table, var_names, basis, z_row):
    headers = ["Базис"] + [v.upper() for v in var_names] + ["B"]
    data = []
    for i in range(table.rows):
        row = [basis[i]] + [str(nsimplify(cell, rational=True)) for cell in table.row(i)]
        data.append(row)
    z = ["Z"] + [str(nsimplify(e, rational=True)) for e in z_row]
    data.append(z)
    print(tabulate(data, headers=headers, tablefmt="grid", stralign="center"))

def choose_pivot_column(c, goal='max'):
    numeric_c = [float(sp.sympify(expr).subs(M, 1000)) for expr in c[:-1]]
    if goal == 'max':
        return numeric_c.index(min(numeric_c)) if min(numeric_c) < 0 else None
    elif goal == 'min':
        pos = [i for i, val in enumerate(numeric_c) if val > 0]
        return max(pos, key=lambda i: numeric_c[i]) if pos else None

def symbolic_simplex_method_final(table, c, goal='max'):
    table = np.array(table.tolist(), dtype=object)
    c = [sp.sympify(val) for val in c]
    num_constraints = table.shape[0]
    num_vars = table.shape[1] - 1
    basis = []

    for col in range(num_vars):
        column = table[:, col]
        if list(column).count(1) == 1 and list(column).count(0) == num_constraints - 1:
            row = list(column).index(1)
            if len(basis) < num_constraints:
                basis.append(col)

    while True:
        pivot_col = choose_pivot_column(c, goal)
        if pivot_col is None:
            break

        ratios = []
        for i in range(num_constraints):
            coeff = table[i, pivot_col]
            rhs = table[i, -1]
            ratio = rhs / coeff if coeff > 0 else float('inf')
            ratios.append(ratio)

        if all(r == float('inf') for r in ratios):
            return None, None, None, None

        pivot_row = ratios.index(min(ratios, key=lambda r: float(r)))
        pivot_element = table[pivot_row, pivot_col]
        table[pivot_row] = [simplify(x / pivot_element) for x in table[pivot_row]]

        for i in range(num_constraints):
            if i != pivot_row:
                table[i] = [simplify(table[i][j] - table[i][pivot_col] * table[pivot_row][j]) for j in range(len(c))]

        c = [simplify(c[j] - c[pivot_col] * table[pivot_row][j]) for j in range(len(c))]
        basis[pivot_row] = pivot_col

    solution = [0] * (len(c) - 1)
    for i in range(num_constraints):
        if basis[i] < len(solution):
            solution[basis[i]] = table[i, -1]
    return table, c, basis, solution

def is_integer(x): return np.isclose(float(x), round(float(x)))
def get_fractional_part(x): return float(x) - np.floor(float(x))

def find_gomory_row_index(basis, solution):
    frac_parts = [get_fractional_part(solution[i]) for i in range(len(solution))]
    var = np.argmax(frac_parts)
    return next(((i, frac_parts[var]) for i, b in enumerate(basis) if b == var), (None, None))

def add_gomory_cut(table, c, solution, basis):
    row_index, _ = find_gomory_row_index(basis, solution)
    if row_index is None: return table, c
    row = table[row_index]
    gomory_row = [-(float(row[i]) - np.floor(float(row[i]))) for i in range(len(row) - 1)]
    gomory_rhs = -(float(row[-1]) - np.floor(float(row[-1])))
    new_row = gomory_row + [1, gomory_rhs]
    new_table = np.hstack((table[:, :-1], np.zeros((table.shape[0], 1)), table[:, -1:]))
    new_table = np.vstack((new_table, np.array(new_row, dtype=object)))
    new_c = np.append(c[:-1], [0, c[-1]])
    basis.append(len(new_c) - 2)
    return new_table, new_c

def dual_simplex(table, c, basis):
    num_constraints, num_vars = table.shape[0], table.shape[1] - 1
    while True:
        if np.all([float(r) >= 0 for r in table[:, -1]]): break
        pivot_row = np.argmin([float(r) for r in table[:, -1]])
        row = table[pivot_row, :-1]
        candidates = [(abs(float(c[j] / row[j])), j) for j in range(num_vars) if float(row[j]) < 0]
        if not candidates: return None, None, None, None
        _, pivot_col = min(candidates)
        table[pivot_row] = [simplify(x / row[pivot_col]) for x in table[pivot_row]]
        for i in range(num_constraints):
            if i != pivot_row:
                table[i] = [simplify(table[i][j] - table[i][pivot_col] * table[pivot_row][j]) for j in range(table.shape[1])]
        c = [simplify(c[j] - c[pivot_col] * table[pivot_row][j]) for j in range(len(c))]
        basis[pivot_row] = pivot_col

    solution = np.zeros(num_vars)
    for i in range(num_constraints):
        if basis[i] < num_vars:
            solution[basis[i]] = table[i, -1]
    return table, c, basis, solution

# === ПРИКЛАД ===
c = [-4, 2]
A = [[1, 1], [-3, 2]]
b = [6, 3]
signs = ["<=", ">="]
goal = 'min'

table, z_row, var_names, basis = build_canonical_table(c, A, b, signs, goal)
print("=== Початкова таблиця ===")
print_table_with_tabulate(table, var_names, basis, z_row)

table, z_row, basis, solution = symbolic_simplex_method_final(table, z_row, goal)
print("\n=== Результат симплекс-методу ===")
print_table_with_tabulate(Matrix(table), var_names, basis, Matrix(z_row))

sol_str = ", ".join([f"{var_names[i]} = {sp.nsimplify(val, rational=True)}" for i, val in enumerate(solution)])
F_value = sum(sp.Rational(c[i]) * solution[i] for i in range(len(c)))
F_value = sp.nsimplify(F_value, rational=True)
print(f"Розв’язок: {sol_str}")
print(f"При цьому значення цільової функції F(x) = {F_value}")

if all(is_integer(val) for val in solution):
    print("\nУсі значення цілі.")
else:
    print("\nПроміжний розв’язок нецілий, додаємо зріз Гоморі.")
    table, z_row = add_gomory_cut(table, z_row, solution, basis)

    print("\n=== Симплекс-таблиця після додавання зрізу Гоморі ===")
    print_table_with_tabulate(Matrix(table), var_names + ["cut"], basis, Matrix(z_row))

    table = np.array([[expr.subs(M, 1000) if isinstance(expr, sp.Basic) else expr for expr in row] for row in table], dtype=float)
    z_row = [expr.subs(M, 1000) if isinstance(expr, sp.Basic) else expr for expr in z_row]

    table, z_row, basis, solution = dual_simplex(table, z_row, basis)

    print("\n=== Після розв’язання методом двоїстого симплексу ===")
    print_table_with_tabulate(Matrix(table), var_names + ["cut"], basis, Matrix(z_row))

    sol_dict = {basis[i]: solution[i] for i in range(len(basis)) if basis[i] < len(var_names)}
    sol_str = ", ".join([f"{var_names[i]} = {sp.nsimplify(sol_dict.get(i, 0), rational=True)}" for i in range(len(var_names))])
    F_value = sum(sp.Rational(c[i]) * sol_dict.get(i, 0) for i in range(len(c)))
    F_value = sp.nsimplify(F_value, rational=True)
    print(f"Новий розв’язок: {sol_str}")
    print(f"F(x) = {F_value}")
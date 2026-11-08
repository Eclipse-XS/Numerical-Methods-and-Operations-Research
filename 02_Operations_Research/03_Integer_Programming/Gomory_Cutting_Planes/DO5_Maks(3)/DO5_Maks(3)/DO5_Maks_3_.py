import sympy as sp
import numpy as np
from sympy import Matrix, nsimplify, simplify
from tabulate import tabulate
from termcolor import colored

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
        basis_name = (
            var_names[basis[i]]
            if isinstance(basis[i], int) and basis[i] < len(var_names)
            else str(basis[i])
        )
        row = [basis_name] + [str(nsimplify(cell, rational=True)) for cell in table.row(i)]
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

def find_gomory_row_index(table):
    fractional_parts = [get_fractional_part(row[-1]) if not is_integer(row[-1]) else 0 for row in table]
    if all(fp == 0 for fp in fractional_parts):
        return None
    return int(np.argmax(fractional_parts))

def add_gomory_cut(table, c, var_names, basis):
    row_index = find_gomory_row_index(table)
    if row_index is None:
        return table, c, var_names, basis

    row = table[row_index]
    gomory_row = [-(get_fractional_part(row[i])) for i in range(len(row) - 1)]
    gomory_rhs = -get_fractional_part(row[-1])

    new_var = f"cut{len(var_names)}"
    gomory_row += [1, gomory_rhs]

    new_table = np.hstack((table[:, :-1], np.zeros((table.shape[0], 1)), table[:, -1:]))
    new_table = np.vstack((new_table, np.array(gomory_row, dtype=object)))
    new_c = np.append(c[:-1], [0, c[-1]])
    var_names.append(new_var)
    basis.append(len(var_names) - 1)
    return new_table, new_c, var_names, basis

def print_highlighted_table(table, c, basis, pivot_row, pivot_col):
    headers = ["Базис"] + [f"x{i+1}" for i in range(table.shape[1] - 1)] + ["B"]
    data = []

    for i, row in enumerate(table):
        row_data = []
        for j, val in enumerate(row):
            val_str = str(nsimplify(val, rational=True))
            if j == pivot_col and i == pivot_row:
                val_str = colored(val_str, "green", attrs=["reverse", "bold"])
            elif j == pivot_col:
                val_str = colored(val_str, "yellow")
            elif i == pivot_row:
                val_str = colored(val_str, "green")
            row_data.append(val_str)
        data.append([f"x{basis[i]+1}"] + row_data)

    z_row = []
    for j, val in enumerate(c):
        val_str = str(nsimplify(val, rational=True))
        if j == pivot_col:
            val_str = colored(val_str, "yellow")
        z_row.append(val_str)
    data.append(["Z"] + z_row)

    print(tabulate(data, headers=headers, tablefmt="grid", stralign="center"))

def dual_simplex(table, c, basis):
    num_constraints, num_vars = table.shape[0], table.shape[1] - 1
    iter_count = 1
    while True:
        if all(table[i, -1] >= 0 for i in range(num_constraints)):
            break

        pivot_row = min(range(num_constraints), key=lambda i: table[i, -1])
        row = table[pivot_row, :-1]
        candidates = [(abs(c[j] / row[j]), j) for j in range(num_vars) if row[j] < 0]
        if not candidates:
            return None, None, None, None

        _, pivot_col = min(candidates, key=lambda x: x[0])
        pivot_element = table[pivot_row, pivot_col]

        print(f"\n=== ІТЕРАЦІЯ {iter_count} ГОМОРІ ===")
        print(f"🔹 Опорний стовпець: {pivot_col + 1}, 🔹 Опорний рядок: {pivot_row + 1}, значення: {nsimplify(pivot_element, rational=True)}")
        print_highlighted_table(table, c, basis, pivot_row, pivot_col)
        iter_count += 1

        table[pivot_row] = [simplify(x / pivot_element) for x in table[pivot_row]]
        for i in range(num_constraints):
            if i != pivot_row:
                factor = table[i, pivot_col]
                table[i] = [simplify(table[i][j] - factor * table[pivot_row][j]) for j in range(table.shape[1])]
        c = [simplify(c[j] - c[pivot_col] * table[pivot_row][j]) for j in range(len(c))]
        basis[pivot_row] = pivot_col

    solution = [0] * num_vars
    for i in range(num_constraints):
        if basis[i] < num_vars:
            solution[basis[i]] = table[i, -1]

    return table, c, basis, solution

# === Приклад ===
c = [1, 3]
A = [
    [8, -5],
    [-6, 5],
    [2, 1],
    [2, 5]
]
b = [40, 60, 14, 10]
signs = ["<=", "<=", "<=", ">="]
goal = 'max'

table, z_row, var_names, basis = build_canonical_table(c, A, b, signs, goal)
print("=== Початкова таблиця ===")
print_table_with_tabulate(table, var_names, basis, z_row)

table, z_row, basis, solution = symbolic_simplex_method_final(table, z_row, goal)
print("\n=== Результат симплекс-методу ===")
print_table_with_tabulate(Matrix(table), var_names, basis, Matrix(z_row))

sol_str = ", ".join([f"{var_names[i]} = {sp.nsimplify(solution[i], rational=True)}" for i in range(len(var_names))])
F_value = sum(sp.Rational(c[i]) * solution[i] for i in range(len(c)))
F_value = sp.nsimplify(F_value, rational=True)
print(f"Розв’язок: {sol_str}")
print(f"F(x) = {F_value}")

if all(is_integer(val) for val in solution):
    print("\nУсі значення цілі.")
else:
    print("\nПроміжний розв’язок нецілий, додаємо зріз Гоморі.")
    table, z_row, var_names, basis = add_gomory_cut(table, z_row, var_names, basis)

    print("\n=== Симплекс-таблиця після додавання зрізу Гоморі ===")
    print_table_with_tabulate(Matrix(table), var_names, basis, Matrix(z_row))

    table, z_row, basis, solution = dual_simplex(table, z_row, basis)

    print("\n=== Після розв’язання методом двоїстого симплексу ===")
    print_table_with_tabulate(Matrix(table), var_names, basis, Matrix(z_row))

    final_solution = [0] * len(var_names)
    for i in range(len(basis)):
        if basis[i] < len(var_names):
            final_solution[basis[i]] = table[i, -1]

    sol_str = ", ".join([f"{var_names[i]} = {sp.nsimplify(final_solution[i], rational=True)}" for i in range(len(var_names))])
    F_value = sum(sp.Rational(c[i]) * final_solution[i] for i in range(len(c)))
    F_value = sp.nsimplify(F_value, rational=True)
    print(f"Новий розв’язок: {sol_str}")
    print(f"F(x) = {F_value}")


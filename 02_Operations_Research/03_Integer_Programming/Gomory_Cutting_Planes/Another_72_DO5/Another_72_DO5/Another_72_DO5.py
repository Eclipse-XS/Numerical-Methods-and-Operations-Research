import sympy as sp
import numpy as np
from sympy import Matrix, nsimplify, simplify, Matrix, Rational, Symbol
from fractions import Fraction

M = sp.Symbol("M")
def is_integer_value(x): return np.isclose(float(x), round(float(x)))
def fractional_part_of(x): return float(x) - np.floor(float(x))

def run_simplex(tableau, cost_row, var_names, base_vars, objective='max'):
    tableau = np.array(tableau.tolist(), dtype=object)
    cost_row = [sp.sympify(val) for val in cost_row]
    num_rows = tableau.shape[0]
    num_cols = tableau.shape[1] - 1  

    while True:
        eval_row = [float(sp.sympify(expr).subs(M, 1000)) for expr in cost_row[:-1]]
        if objective == 'max':
            if min(eval_row) >= 0:
                break
            pivot_col = eval_row.index(min(eval_row))
        elif objective == 'min':
            candidates = [i for i, v in enumerate(eval_row) if v > 0]
            if not candidates:
                break
            pivot_col = max(candidates, key=lambda i: eval_row[i])
        ratios = []
        for i in range(num_rows):
            coefficient = tableau[i, pivot_col]
            rhs = tableau[i, -1]
            ratio = rhs / coefficient if coefficient > 0 else float('inf')
            ratios.append(ratio)

        if all(r == float('inf') for r in ratios):
            return None, None, None, None

        pivot_row = ratios.index(min(ratios, key=float))
        pivot_elem = tableau[pivot_row, pivot_col]
        tableau[pivot_row] = [simplify(x / pivot_elem) for x in tableau[pivot_row]]
        for i in range(num_rows):
            if i != pivot_row:
                factor = tableau[i, pivot_col]
                tableau[i] = [simplify(tableau[i][j] - factor * tableau[pivot_row][j]) for j in range(tableau.shape[1])]

        factor = cost_row[pivot_col]
        cost_row = [simplify(cost_row[j] - factor * tableau[pivot_row][j]) for j in range(len(cost_row))]
        base_vars[pivot_row] = var_names[pivot_col]

    solution = [0] * num_cols
    for i in range(num_rows):
        if base_vars[i] in var_names[:num_cols]:
            idx = var_names.index(base_vars[i])
            solution[idx] = tableau[i, -1]

    return tableau, cost_row, base_vars, solution

def cut(tableau, cost_row, variable_names, base_vars):
    target_vars = ["x1", "x2"]
    row_candidates = []

    for row_idx, base_var in enumerate(base_vars):
        if base_var in target_vars:
            value = tableau[row_idx, -1]
            if not is_integer_value(value):
                frac_part = fractional_part_of(value)
                row_candidates.append((row_idx, frac_part))

    if not row_candidates:
        return tableau, cost_row, variable_names, base_vars, False

    selected_row, max_fraction = max(row_candidates, key=lambda x: x[1])
    if max_fraction == 0:
        return tableau, cost_row, variable_names, base_vars, False

    row = tableau[selected_row]
    cut_row = [Fraction(-fractional_part_of(row[i])).limit_denominator(1000) for i in range(len(row) - 1)]
    cut_row += [1, Fraction(-fractional_part_of(row[-1])).limit_denominator(1000)]

    extended_tableau = np.hstack((tableau[:, :-1], np.zeros((tableau.shape[0], 1)), tableau[:, -1:]))
    extended_tableau = np.vstack((extended_tableau, np.array(cut_row, dtype=object)))

    extended_cost_row = np.append(cost_row[:-1], [0, cost_row[-1]])
    new_var = f"cut{sum(1 for v in variable_names if v.startswith('cut')) + 1}"
    variable_names.append(new_var)
    base_vars.append(new_var)

    return extended_tableau, extended_cost_row, variable_names, base_vars, True

def run_dual(table, c, var_names, basis):
    num_constraints, num_vars = table.shape[0], table.shape[1] - 1
    while True:
        if all(table[i, -1] >= 0 for i in range(num_constraints)):
            break
        pivot_row = min((i for i in range(num_constraints) if table[i, -1] < 0), key=lambda i: table[i, -1])
        row = table[pivot_row, :-1]
        candidates = [(abs(c[j] / row[j]), j) for j in range(num_vars) if row[j] < 0]
        if not candidates:
            return None, None, None, None
        _, pivot_col = min(candidates, key=lambda x: float(x[0].subs(M, 1000) if isinstance(x[0], sp.Basic) else x[0]))
        pivot_element = table[pivot_row, pivot_col]
        table[pivot_row] = [simplify(x / pivot_element) for x in table[pivot_row]]
        for i in range(num_constraints):
            if i != pivot_row:
                factor = table[i, pivot_col]
                table[i] = [simplify(table[i][j] - factor * table[pivot_row][j]) for j in range(table.shape[1])]
        factor = c[pivot_col]
        c = [simplify(c[j] - factor * table[pivot_row][j]) for j in range(len(c))]
        basis[pivot_row] = var_names[pivot_col]
    solution = [0] * num_vars
    for i in range(num_constraints):
        if basis[i] in var_names[:num_vars]:
            idx = var_names.index(basis[i])
            solution[idx] = table[i, -1]
    return table, c, basis, solution


def generate_canonical_table(obj_coeffs, constraints, rhs_values, constraint_types, objective='max'):
    num_vars = len(obj_coeffs)
    A, b= Matrix(constraints), Matrix(rhs_values).reshape(len(rhs_values), 1)
    tableau = A
    var_names = [f"x{i+1}" for i in range(num_vars)]
    basis = []
    artificial_indices = []

    next_idx = num_vars + 1

    for i, sign in enumerate(constraint_types):
        if sign in ("<=", ">="):
            slack = Matrix.zeros(len(rhs_values), 1)
            slack[i] = 1 if sign == "<=" else -1
            tableau = tableau.row_join(slack)
            var_names.append(f"x{next_idx}")
            if sign == "<=":
                basis.append(f"x{next_idx}")
            next_idx += 1

        if sign in (">=", "="):
            artificial = Matrix.zeros(len(rhs_values), 1)
            artificial[i] = 1
            tableau = tableau.row_join(artificial)
            w_name = f"w{len([v for v in var_names if v.startswith('w')]) + 1}"
            var_names.append(w_name)
            basis.append(w_name)
            artificial_indices.append((i, tableau.cols - 1))

    tableau = tableau.row_join(b)

    cost_row = Matrix([0] * (tableau.cols - 1) + [0])
    for i, val in enumerate(obj_coeffs):
        cost_row[i] = Rational(val)

    mult = 1 if objective == 'min' else -1
    for row_i, _ in artificial_indices:
        row = tableau.row(row_i)
        for j in range(len(cost_row) - 1):
            cost_row[j] += mult * M * (-row[j])
        cost_row[-1] += mult * M * row[-1]

    for _, col in artificial_indices:
        cost_row[col] = 0

    for j in range(len(cost_row) - 1):
        cost_row[j] *= -1

    return tableau, cost_row, var_names, basis

A = [ [1, 1], [-5, 1], [1, -5], [1, 1]]
b = [1, 0, 0, 6]
c = [2, 2]
signs = ['>=', '<=', '<=', '<=']
goal = 'min'

table, z_row, var_names, basis = generate_canonical_table(c, A, b, signs, goal)
table, z_row, basis, solution = run_simplex(table, z_row, var_names, basis, goal)

solution_str = ", ".join([f"{var_names[i]} = {nsimplify(solution[i], rational=True)}"
                          for i in range(len(var_names)) if var_names[i] in ["x1", "x2"]])
F_value = sum(sp.Rational(c[i]) * solution[i] for i in range(len(c)))
F_value = sp.nsimplify(F_value, rational=True)

print("═" * 50)
print("Результат обчислення симплекс-методу:")
print(f"    ({solution_str})")
print(f"    F(x) = {F_value}")
print("═" * 50)

cut_count = 0
if not all(is_integer_value(solution[i]) for i in range(len(c))):
    print("Розв'язок не є цілочисельний, проводимо зріз Гоморі")

while not all(is_integer_value(solution[i]) for i in range(len(c))):
    cut_count += 1
    print("\n" + "-" * 50)
    print(f"Ітерація {cut_count} — Додавання зрізу Гоморі")
    print("-" * 50)

    table, z_row, var_names, basis, added = cut(table, z_row, var_names, basis)
    if not added:
        print("Неможливо додати новий зріз.")
        break

    table, z_row, basis, solution = run_dual(table, z_row, var_names, basis)

    solution_str = ", ".join([f"{var_names[i]} = {nsimplify(solution[i], rational=True)}"
                              for i in range(len(var_names)) if var_names[i] in ["x1", "x2"]])
    F_value = sum(sp.Rational(c[i]) * solution[i] for i in range(len(c)))
    F_value = sp.nsimplify(F_value, rational=True)

    print("Результат обчислення двоїстого симплекс-методу:")
    print(f"    ({solution_str})")
    print(f"    F(x) = {F_value}")

    if all(is_integer_value(solution[i]) for i in range(len(c))):
        print("Результат є: цілочисельний — зупиняємось.")
    else:
        print("Результат є: нецілочисельний — проводимо ще один зріз.")

print("\n" + "═" * 50)
print("Фінальний розв’язок:")
sol_str = ", ".join([f"{var_names[i]} = {nsimplify(solution[i], rational=True)}"
                     for i in range(len(var_names)) if var_names[i] in ["x1", "x2"]])
print(f"    ({sol_str})")
print(f"    F(x) = {F_value}")
print(f"    Кількість зрізів Гоморі: {cut_count}")
print("═" * 50)
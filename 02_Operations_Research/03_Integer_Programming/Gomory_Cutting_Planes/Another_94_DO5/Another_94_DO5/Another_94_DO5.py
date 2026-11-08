import sympy as sp
import numpy as np
from sympy import Matrix, nsimplify, simplify

M = sp.Symbol("M")

objective_coeffs = [1, 1]
constraint_matrix = [
    [1, 1],
    [3, 7],
    [1, 2],
    [1, 0],
    [0, 1],
]
right_hand_side = [8, 21, 6, 7, 7]
constraint_signs = ["<=", ">=", ">=", "<=", "<="]
optimization_type = "min"

def setup_initial_table(obj_coeffs, constraints_matrix, rhs_vector, constraint_signs, optimization_goal='max'):
    num_variables = len(obj_coeffs)
    constraints_matrix = Matrix(constraints_matrix)
    rhs_vector = Matrix(rhs_vector).reshape(len(rhs_vector), 1)
    variable_labels = [f"x{i+1}" for i in range(num_variables)]
    basic_variables = []
    artificial_entries = []

    tableau = constraints_matrix
    next_var_idx = num_variables + 1

    for idx, sign in enumerate(constraint_signs):
        if sign == "<=" :
            slack = Matrix.zeros(len(constraint_signs), 1)
            slack[idx] = 1
            tableau = tableau.row_join(slack)
            slack_var = f"x{next_var_idx}"
            next_var_idx += 1
            variable_labels.append(slack_var)
            basic_variables.append(slack_var)

        elif sign == ">=":
            surplus = Matrix.zeros(len(constraint_signs), 1)
            surplus[idx] = -1
            tableau = tableau.row_join(surplus)
            surplus_var = f"x{next_var_idx}"
            next_var_idx += 1
            variable_labels.append(surplus_var)
            artificial = Matrix.zeros(len(constraint_signs), 1)
            artificial[idx] = 1
            tableau = tableau.row_join(artificial)
            art_var = f"a{len([v for v in variable_labels if v.startswith('a')]) + 1}"
            variable_labels.append(art_var)
            basic_variables.append(art_var)
            artificial_entries.append((idx, tableau.cols - 1))

        elif sign == "=":
            artificial = Matrix.zeros(len(constraint_signs), 1)
            artificial[idx] = 1
            tableau = tableau.row_join(artificial)
            art_var = f"a{len([v for v in variable_labels if v.startswith('a')]) + 1}"
            variable_labels.append(art_var)
            basic_variables.append(art_var)
            artificial_entries.append((idx, tableau.cols - 1))

    tableau = tableau.row_join(rhs_vector)

    objective_row = Matrix([0] * (tableau.cols - 1) + [0])
    for i in range(len(obj_coeffs)):
        objective_row[i] = sp.Rational(obj_coeffs[i])

    multiplier = 1 if optimization_goal == 'min' else -1
    for row_idx, _ in artificial_entries:
        row = tableau.row(row_idx)
        rhs = row[-1]
        for j in range(len(objective_row) - 1):
            objective_row[j] += multiplier * M * (-row[j])
        objective_row[-1] += multiplier * M * rhs

    for _, art_col in artificial_entries:
        objective_row[art_col] = 0

    for i in range(len(objective_row) - 1):
        objective_row[i] *= -1

    return tableau, objective_row, variable_labels, basic_variables

def run_simplex(table, c, var_names, basis, goal='max'):
    num_vars = c.shape[0]
    num_constraints = table.shape[0]

    def select_pivot(c):
        numeric_c = [float(expr.subs(M, 1000)) for expr in c[:-1]]
        if goal == 'max':
            return numeric_c.index(min(numeric_c)) if min(numeric_c) < 0 else None
        elif goal == 'min':
            pos = [i for i, val in enumerate(numeric_c) if val > 0]
            return max(pos, key=lambda i: numeric_c[i]) if pos else None

    while True:
        pivot_col = select_pivot(c)
        if pivot_col is None:
            break

        ratios = []
        for i in range(num_constraints):
            coeff = table[i, pivot_col]
            rhs = table[i, -1]
            ratio = rhs / coeff if coeff > 0 else float('inf')
            ratios.append(ratio)

        if all(r == float('inf') for r in ratios):
            return None, None, None

        pivot_row = ratios.index(min(ratios, key=lambda r: float(r)))
        pivot_element = table[pivot_row, pivot_col]

        table.row_op(pivot_row, lambda x, _: simplify(x / pivot_element))

        for i in range(num_constraints):
            if i != pivot_row:
                factor = table[i, pivot_col]
                table.row_op(i, lambda x, j: simplify(x - factor * table[pivot_row, j]))

        factor = c[pivot_col]
        c = Matrix([simplify(c[j] - factor * table[pivot_row, j]) for j in range(num_vars)])
        basis[pivot_row] = var_names[pivot_col]

    solution = [0] * num_vars
    for i in range(num_constraints):
        col = var_names.index(basis[i])
        solution[col] = table[i, -1]

    return solution, c[-1], basis

def is_integer(x):
    return np.isclose(float(x), round(float(x)))

def fractional_part(x):
    return float(x) - np.floor(float(x))

def find_gomory_cutting_row(table):
    fractions = [fractional_part(row[-1]) if not is_integer(row[-1]) else 0 for row in table]
    return None if all(fp == 0 for fp in fractions) else int(np.argmax(fractions))

def add_gomory_cut(table, cost_vector, var_names, basis):
    target_row = find_gomory_cutting_row(table)
    if target_row is None:
        return table, cost_vector, var_names, basis, ""
    row = table[target_row]
    cut_vector = [-fractional_part(row[i]) for i in range(len(row) - 1)] + [1, -fractional_part(row[-1])]
    new_var = f"x{len(var_names)+1}"
    readable_cut = " + ".join(
        [f"{coef:+.2f}*{var_names[i]}" for i, coef in enumerate(cut_vector[:-2]) if abs(coef) > 1e-8]
    ) + f" + {cut_vector[-2]:+.2f}*{new_var} = {cut_vector[-1]:.2f}"
    expanded_table = np.hstack((table[:, :-1], np.zeros((table.shape[0], 1)), table[:, -1:]))
    expanded_table = np.vstack((expanded_table, np.array(cut_vector, dtype=object)))
    updated_cost = np.append(cost_vector[:-1], [0, cost_vector[-1]])
    var_names.append(new_var)
    basis.append(len(var_names) - 1)
    return expanded_table, updated_cost, var_names, basis, readable_cut

def run_dual_simplex_method(table, cost_vector, basis):
    num_rows, num_cols = table.shape[0], table.shape[1] - 1
    while any(table[i, -1] < 0 for i in range(num_rows)):
        pivot_row = min(range(num_rows), key=lambda i: table[i, -1])
        row = table[pivot_row, :-1]
        candidates = [(abs(cost_vector[j] / row[j]), j) for j in range(num_cols) if row[j] < 0]
        if not candidates:
            return None, None, None, None
        _, pivot_col = min(candidates, key=lambda x: x[0])
        pivot_elem = table[pivot_row, pivot_col]
        table[pivot_row] = [simplify(x / pivot_elem) for x in table[pivot_row]]
        for i in range(num_rows):
            if i != pivot_row:
                factor = table[i, pivot_col]
                table[i] = [simplify(table[i][j] - factor * table[pivot_row][j]) for j in range(table.shape[1])]
        cost_vector = [
            simplify(cost_vector[j] - cost_vector[pivot_col] * table[pivot_row][j])
            for j in range(len(cost_vector))
        ]
        basis[pivot_row] = pivot_col

    solution = [0] * num_cols
    for i in range(num_rows):
        if basis[i] < num_cols:
            solution[basis[i]] = table[i, -1]
    return table, cost_vector, basis, solution

print(f"› Цільова функція: F(x1, x2) = {objective_coeffs[0]}*x1 + {objective_coeffs[1]}*x2")
print(f"› Тип оптимізації: {optimization_type.upper()}\n")
print("› Обмеження:")
for i, row in enumerate(constraint_matrix):
    print(f"  {row[0]}*x1 + {row[1]}*x2 {constraint_signs[i]} {right_hand_side[i]}")

simplex_table, cost_row, all_vars, current_basis = setup_initial_table(
    objective_coeffs, constraint_matrix, right_hand_side, constraint_signs, optimization_type
)
preliminary_solution, z_value, current_basis = run_simplex(
    simplex_table.copy(), cost_row.copy(), all_vars.copy(), current_basis.copy(), optimization_type
)

print("\n--- Етап I: Симплекс-метод ---")
print(f"› ( x1 = {nsimplify(preliminary_solution[0], rational=True)} ; x2 = {nsimplify(preliminary_solution[1], rational=True)} )")
print(f"› Значення F(x1;x2) = {nsimplify(z_value, rational=True)}")

if all(is_integer(val) for val in preliminary_solution):
    print("\nРозв’язок є цілочисельним. Метод Гоморі не застосовується.")
else:
    print("\nРозв’язок не є цілочисельним → застосовується метод Гоморі...")
    simplex_table, cost_row, all_vars, current_basis, added_cut = add_gomory_cut(
        simplex_table, cost_row, all_vars, current_basis
    )
    print(f"› Додане обмеження:\n  {added_cut}")
    simplex_table, cost_row, current_basis, final_solution = run_dual_simplex_method(
        simplex_table, cost_row, current_basis
    )
    final_objective = sum(sp.Rational(objective_coeffs[i]) * final_solution[i] for i in range(len(objective_coeffs)))
    print("\n--- Етап II: Двоїстий симплекс-метод ---")
    print(f"› ( x1 = {nsimplify(final_solution[0], rational=True)} ; x2 = {nsimplify(final_solution[1], rational=True)} )")
    print(f"› Підсумкове значення F(x1; x2) = {nsimplify(final_objective, rational=True)}")
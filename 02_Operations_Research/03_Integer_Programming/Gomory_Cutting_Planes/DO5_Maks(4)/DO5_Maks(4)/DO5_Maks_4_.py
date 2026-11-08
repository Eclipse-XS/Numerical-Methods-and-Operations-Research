import sympy as sp, numpy as np
from sympy import Matrix, nsimplify, simplify
M = sp.Symbol("M")
objective_coeffs = [1, 3]
constraint_matrix = [[8, -5], [-6, 5], [2, 1], [2, 5]]
right_hand_side = [40, 60, 14, 10]
constraint_signs = ["<=", "<=", "<=", ">="]
optimization_type = 'min'

def create_initial_table(c, A, b, signs, target='max'):
    A, b = Matrix(A), Matrix(b).reshape(len(b), 1)
    table = A.row_join(Matrix.eye(len(b)) * 0)
    variable_names = [f"x{i+1}" for i in range(len(c))]
    slack_counter = 0
    artificial_counter = 0
    basis_variables = []
    artificial_entries = []
    for i, sign in enumerate(signs):
        if sign in ["<=", ">="]:
            multiplier = 1 if sign == "<=" else -1
            table[i, len(c) + slack_counter] = multiplier
            variable_names.append(f"x{len(variable_names)+1}")
            if sign == "<=" :
                basis_variables.append(variable_names[-1])
            else:
                table = table.row_join(Matrix.zeros(table.rows, 1))
                table[i, table.cols - 1] = 1
                variable_names.append(f"x{len(variable_names)+1}")
                basis_variables.append(variable_names[-1])
                artificial_entries.append((i, variable_names[-1]))
                artificial_counter += 1
            slack_counter += 1
        elif sign == "=":
            table = table.row_join(Matrix.zeros(table.rows, 1))
            table[i, table.cols - 1] = 1
            variable_names.append(f"x{len(variable_names)+1}")
            basis_variables.append(variable_names[-1])
            artificial_entries.append((i, variable_names[-1]))
            artificial_counter += 1
    table = table.row_join(b)
    cost_row = [sp.Rational(c[i]) for i in range(len(c))] + [0] * (table.cols - len(c) - 1) + [0]
    sign_mult = 1 if target == 'min' else -1
    for row_idx, art_var in artificial_entries:
        rhs = table[row_idx, -1]
        for j in range(table.cols - 1):
            if variable_names[j] != art_var:
                cost_row[j] += sign_mult * M * (-table[row_idx, j])
        cost_row[variable_names.index(art_var)] = 0
        cost_row[-1] += sign_mult * M * rhs
    return table, Matrix([-val for val in cost_row]), variable_names, basis_variables

def select_pivot_column(cost_vector, target='max'):
    evaluations = [float(sp.sympify(expr).subs(M, 1000)) for expr in cost_vector[:-1]]
    return (evaluations.index(min(evaluations)) if min(evaluations) < 0 else None) if target == 'max' else (max((i for i, val in enumerate(evaluations) if val > 0), default=None, key=lambda i: evaluations[i]))

def perform_simplex(table, cost_vector, target='max'):
    table, cost_vector = np.array(table.tolist(), dtype=object), [sp.sympify(val) for val in cost_vector]
    num_rows, num_cols = table.shape[0], table.shape[1] - 1
    basis = [col for col in range(num_cols) if list(table[:, col]).count(1) == 1 and list(table[:, col]).count(0) == num_rows - 1][:num_rows]
    while (pivot_col := select_pivot_column(cost_vector, target)) is not None:
        ratios = [table[i, -1] / table[i, pivot_col] if table[i, pivot_col] > 0 else float('inf') for i in range(num_rows)]
        if all(r == float('inf') for r in ratios): return None, None, None, None
        pivot_row = ratios.index(min(ratios, key=float))
        pivot_elem = table[pivot_row, pivot_col]
        table[pivot_row] = [simplify(x / pivot_elem) for x in table[pivot_row]]
        for i in range(num_rows):
            if i != pivot_row:
                table[i] = [simplify(table[i][j] - table[i][pivot_col] * table[pivot_row][j]) for j in range(len(cost_vector))]
        cost_vector = [simplify(cost_vector[j] - cost_vector[pivot_col] * table[pivot_row][j]) for j in range(len(cost_vector))]
        basis[pivot_row] = pivot_col
    solution = [0] * (len(cost_vector) - 1)
    for i in range(num_rows):
        if basis[i] < len(solution):
            solution[basis[i]] = table[i, -1]
    return table, cost_vector, basis, solution

def is_whole_number(x): return np.isclose(float(x), round(float(x)))
def get_decimal_part(x): return float(x) - np.floor(float(x))
def locate_cutting_row(table):
    fractions = [get_decimal_part(row[-1]) if not is_whole_number(row[-1]) else 0 for row in table]
    return None if all(fp == 0 for fp in fractions) else int(np.argmax(fractions))

def apply_gomory_cut(table, cost_vector, var_names, basis):
    target_row = locate_cutting_row(table)
    if target_row is None: return table, cost_vector, var_names, basis, ""
    row = table[target_row]
    cut_vector = [-get_decimal_part(row[i]) for i in range(len(row) - 1)] + [1, -get_decimal_part(row[-1])]
    new_var = f"x{len(var_names)+1}"
    readable_cut = " + ".join([f"{coef:+.2f}*{var_names[i]}" for i, coef in enumerate(cut_vector[:-2]) if abs(coef) > 1e-8]) + f" + {cut_vector[-2]:+.2f}*{new_var} = {cut_vector[-1]:.2f}"
    expanded_table = np.hstack((table[:, :-1], np.zeros((table.shape[0], 1)), table[:, -1:]))
    expanded_table = np.vstack((expanded_table, np.array(cut_vector, dtype=object)))
    updated_cost = np.append(cost_vector[:-1], [0, cost_vector[-1]])
    var_names.append(new_var)
    basis.append(len(var_names) - 1)
    return expanded_table, updated_cost, var_names, basis, readable_cut

def apply_dual_simplex(table, cost_vector, basis):
    num_rows, num_cols = table.shape[0], table.shape[1] - 1
    while any(table[i, -1] < 0 for i in range(num_rows)):
        pivot_row = min(range(num_rows), key=lambda i: table[i, -1])
        row = table[pivot_row, :-1]
        candidates = [(abs(cost_vector[j] / row[j]), j) for j in range(num_cols) if row[j] < 0]
        if not candidates: return None, None, None, None
        _, pivot_col = min(candidates, key=lambda x: x[0])
        pivot_elem = table[pivot_row, pivot_col]
        table[pivot_row] = [simplify(x / pivot_elem) for x in table[pivot_row]]
        for i in range(num_rows):
            if i != pivot_row:
                factor = table[i, pivot_col]
                table[i] = [simplify(table[i][j] - factor * table[pivot_row][j]) for j in range(table.shape[1])]
        cost_vector = [simplify(cost_vector[j] - cost_vector[pivot_col] * table[pivot_row][j]) for j in range(len(cost_vector))]
        basis[pivot_row] = pivot_col
    solution = [0] * num_cols
    for i in range(num_rows):
        if basis[i] < num_cols:
            solution[basis[i]] = table[i, -1]
    return table, cost_vector, basis, solution

print(f"F = {objective_coeffs[0]}*x1 + {objective_coeffs[1]}*x2")
print(f"Обраний тип оптимізації: {optimization_type.upper()}\n")
print("Система умов:")
for i, row in enumerate(constraint_matrix):
    print(f"{row[0]}*x1 + {row[1]}*x2 {constraint_signs[i]} {right_hand_side[i]}")

simplex_table, cost_row, all_vars, current_basis = create_initial_table(objective_coeffs, constraint_matrix, right_hand_side, constraint_signs, optimization_type)
simplex_table, cost_row, current_basis, preliminary_solution = perform_simplex(simplex_table, cost_row, optimization_type)

print("\nОбчислення симплекс-методом:")
print(f"x1 = {nsimplify(preliminary_solution[0], rational=True)}, x2 = {nsimplify(preliminary_solution[1], rational=True)}")
objective_value = sum(sp.Rational(objective_coeffs[i]) * preliminary_solution[i] for i in range(len(objective_coeffs)))
print(f"F = {nsimplify(objective_value, rational=True)}")

if all(is_whole_number(val) for val in preliminary_solution):
    print("\nРозв’язок є цілочисельним. Метод Гоморі не застосовується.")
else:
    print("\nОскільки (x1; x2) не цілочисельний, застосовується метод Гоморі.")
    simplex_table, cost_row, all_vars, current_basis, added_cut = apply_gomory_cut(simplex_table, cost_row, all_vars, current_basis)
    print(f"Додано обмеження: \n{added_cut}")
    simplex_table, cost_row, current_basis, final_solution = apply_dual_simplex(simplex_table, cost_row, current_basis)
    final_objective = sum(sp.Rational(objective_coeffs[i]) * final_solution[i] for i in range(len(objective_coeffs)))
    print("\nРезультат після зрізу Гоморі:")
    print(f"x1 = {nsimplify(final_solution[0], rational=True)}, x2 = {nsimplify(final_solution[1], rational=True)}")
    print(f"F = {nsimplify(final_objective, rational=True)}")
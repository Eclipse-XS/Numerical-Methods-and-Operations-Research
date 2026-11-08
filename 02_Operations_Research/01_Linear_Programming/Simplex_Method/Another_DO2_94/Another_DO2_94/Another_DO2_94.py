import sympy as sp
from sympy import Matrix, simplify
M = sp.Symbol("M")

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

def solve_and_report(A, b, c, signs, goal):
    table, z_row, var_names, basis_vars = setup_initial_table(c, A, b, signs, goal)
    solution, optimal_value, final_basis = run_simplex(table.copy(), z_row.copy(), var_names, basis_vars.copy(), goal)
    print(f"\nРезультат {goal}")
    for i, val in enumerate(solution):
        if val != 0:
            val = float(val)
            val = round(val, 6)
            print(f"{var_names[i]} = {val}")
    optimal_value = float(optimal_value)
    optimal_value = round(optimal_value, 6)
    print(f"Оптимальне значення цільової функції: f(x) = {optimal_value}")

A = [
    [1, 1],
    [3, 7],
    [1, 2],
    [1, 0],
    [0, 1]
]
b = [8, 21, 6, 7.5, 7.5]
c = [1, 1]
signs = ['<=', '>=', '>=', '<=', '<=']

solve_and_report(A, b, c, signs, 'max')
solve_and_report(A, b, c, signs, 'min')
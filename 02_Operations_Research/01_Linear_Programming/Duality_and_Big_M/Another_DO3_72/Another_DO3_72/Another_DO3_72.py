from tabulate import tabulate
import numpy as np
import sympy as sp
from sympy import Matrix, simplify

M = sp.Symbol("M")

def preprocess_constraints(matrix, rhs, signs, objective_type, invert=False):
    target_sign = '>=' if objective_type == 'min' else '<='
    new_matrix, new_rhs, new_signs = [], [], []

    for i in range(len(rhs)):
        row, value, sign = matrix[i], rhs[i], signs[i]
        if (invert and value < 0) or (not invert and sign != target_sign):
            row = [-x for x in row]
            value *= -1
            sign = {'>=': '<=', '<=': '>=', '=': '='}[sign]
        new_matrix.append(row)
        new_rhs.append(value)
        new_signs.append(sign)

    return new_matrix, new_rhs, new_signs

def construct_dual(matrix, rhs, cost, direction='min'):
    matrix = np.array(matrix)
    rhs = np.array(rhs)
    cost = np.array(cost)

    return {
        'type': 'max' if direction == 'min' else 'min',
        'A': matrix.T.tolist(),
        'b': cost.tolist(),
        'c': rhs.tolist(),
        'signs': ['<=' if direction == 'min' else '>='] * len(matrix[0])
    }

def format_expression(coefficients, var='x'):
    return " ".join([
        f"{'' if coef < 0 or i == 0 else '+'}{coef}{var}{i + 1}"
        for i, coef in enumerate(coefficients)
    ])

def display_primal(c, A, b, goal):
    print(f"\n{'MIN(F)' if goal == 'min' else 'MAX(F)'} = {format_expression(c)}")
    print("Subject to:")
    constraints = [
        [format_expression(A[i]), '>=', b[i]]
        for i in range(len(A))
    ]
    print(tabulate(constraints, headers=["Constraint", "", "Value"], tablefmt="fancy_grid"))

def display_dual(dual):
    print(f"\n{'MAX(Z)' if dual['type'] == 'max' else 'MIN(Z)'} = {format_expression(dual['c'], 'y')}")
    print("Subject to:")
    constraints = [
        [format_expression(dual['A'][i], 'y'), dual['signs'][i], dual['b'][i]]
        for i in range(len(dual['A']))
    ]
    print(tabulate(constraints, headers=["Constraint", "", "Value"], tablefmt="fancy_grid"))

def display_solution(sol, opt):
    print("\nSOLUTION:")
    results = [
        [f"y{i+1}", simplify(val)]
        for i, val in enumerate(sol) if val != 0
    ]
    results.append(["Z", simplify(opt)])
    print(tabulate(results, headers=["Variable", "Value"], tablefmt="fancy_grid"))

def create_simplex_table(c, A, b, signs, goal='min'):
    A = Matrix(A)
    b = Matrix(b).reshape(len(b), 1)
    variable_names = [f"y{i+1}" for i in range(len(c))]
    tableau = A
    basis = []
    artificial_vars = []

    for i, sign in enumerate(signs):
        if sign in ['<=', '≥']:
            slack = Matrix.zeros(len(b), 1)
            slack[i] = 1 if sign == '<=' else -1
            tableau = tableau.row_join(slack)
            variable_names.append(f"s{i+1}")
            if sign == '<=':
                basis.append(variable_names[-1])
        if sign in ['>=', '=', '≥']:
            art = Matrix.zeros(len(b), 1)
            art[i] = 1
            tableau = tableau.row_join(art)
            variable_names.append(f"a{i+1}")
            basis.append(variable_names[-1])
            artificial_vars.append((i, tableau.shape[1] - 1))

    tableau = tableau.row_join(b)
    z_row = Matrix([0] * (tableau.shape[1] - 1) + [0])

    for i in range(len(c)):
        z_row[i] = sp.Rational(c[i])

    for i, _ in artificial_vars:
        row = tableau.row(i)
        rhs = tableau[i, -1]
        for j in range(len(z_row) - 1):
            z_row[j] += (M if goal == 'min' else -M) * (-row[j])
        z_row[-1] += (M if goal == 'min' else -M) * rhs

    for i in range(len(z_row) - 1):
        z_row[i] *= -1
    for _, col in artificial_vars:
        z_row[col] = 0

    return tableau, z_row, variable_names, basis

def simplex_algorithm(table, z, names, basis, goal='min'):
    while True:
        z_numeric = [float(expr.subs(M, 1000)) for expr in z[:-1]]
        pivot_col = (max if goal == 'min' else min)(
            (i for i, val in enumerate(z_numeric) if (val > 0 if goal == 'min' else val < 0)),
            default=None,
            key=lambda i: z_numeric[i]
        )
        if pivot_col is None:
            break

        ratios = [table[i, -1] / table[i, pivot_col] if table[i, pivot_col] > 0 else float('inf') for i in range(table.shape[0])]
        if all(r == float('inf') for r in ratios):
            return None, None, None

        pivot_row = ratios.index(min(ratios))
        pivot_value = table[pivot_row, pivot_col]
        table.row_op(pivot_row, lambda val, _: simplify(val / pivot_value))

        for i in range(table.shape[0]):
            if i != pivot_row:
                factor = table[i, pivot_col]
                table.row_op(i, lambda val, j: simplify(val - factor * table[pivot_row, j]))

        factor = z[pivot_col]
        z = Matrix([simplify(z[j] - factor * table[pivot_row, j]) for j in range(z.shape[0])])
        basis[pivot_row] = names[pivot_col]

    solution = [0] * z.shape[0]
    for i, base_var in enumerate(basis):
        if base_var in names:
            solution[names.index(base_var)] = table[i, -1]

    return solution, z[-1], basis

# === INPUT DATA ===
constraints_matrix = [ [1, 1], [-5, 1], [-1, 5], [1, 1]]
rhs_vector = [1, 10, 10, 6]
objective_coeffs = [2, 2]
constraint_signs = ['>=', '<=', '>=', '<=']
optimization_type = 'max'
A_adj, b_adj, signs_adj = preprocess_constraints(constraints_matrix, rhs_vector, constraint_signs, optimization_type)
display_primal(objective_coeffs, A_adj, b_adj, optimization_type)
dual = construct_dual(A_adj, b_adj, objective_coeffs, optimization_type)
display_dual(dual)
dual['A'], dual['b'], dual['signs'] = preprocess_constraints(dual['A'], dual['b'], dual['signs'], dual['type'], invert=True)
tableau, z_row, var_names, base_vars = create_simplex_table(dual['c'], dual['A'], dual['b'], dual['signs'], dual['type'])
solution, optimal, _ = simplex_algorithm(tableau, z_row, var_names, base_vars, dual['type'])
if solution:
    display_solution(solution, optimal)
else:
    print("No solution found or problem is unbounded.")
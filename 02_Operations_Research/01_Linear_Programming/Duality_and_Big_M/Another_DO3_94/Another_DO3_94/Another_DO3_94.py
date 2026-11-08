import numpy as np
import sympy as sp
from sympy import Matrix, simplify
from tabulate import tabulate

M = sp.Symbol("M")

def normalize_primal(A, b, signs, sense):
    A_new, b_new, signs_new = [], [], []
    target_sign = '>=' if sense == 'min' else '<='
    for row, rhs, sign in zip(A, b, signs):
        if sign != target_sign:
            row = [-a for a in row]
            rhs = -rhs
            sign = target_sign
        A_new.append(row)
        b_new.append(rhs)
        signs_new.append(sign)
    return A_new, b_new, signs_new

def fix_negative_rhs_for_dual(A, b, signs):
    A_fixed, b_fixed, signs_fixed = [], [], []
    for i in range(len(b)):
        if b[i] < 0:
            A_fixed.append([-x for x in A[i]])
            b_fixed.append(-b[i])
            new_sign = '<=' if signs[i] == '>=' else '>=' if signs[i] == '<=' else '='
            signs_fixed.append(new_sign)
        else:
            A_fixed.append(A[i])
            b_fixed.append(b[i])
            signs_fixed.append(signs[i])
    return A_fixed, b_fixed, signs_fixed

def build_dual(A, b, c, sense='min'):
    A, b, c = np.array(A), np.array(b), np.array(c)
    A_dual = A.T.tolist()
    c_dual = b.tolist()
    b_dual = c.tolist()
    dual_sense = 'max' if sense == 'min' else 'min'
    constraint_signs = ['<=' if sense == 'min' else '>='] * len(A_dual)
    return {
        'dual_sense': dual_sense,
        'A_dual': A_dual,
        'b_dual': b_dual,
        'c_dual': c_dual,
        'constraint_signs': constraint_signs
    }

def format_expression(coeffs, vars_prefix='x'):
    return " ".join([
        f"{'' if coef < 0 or i == 0 else '+'}{coef}{vars_prefix}{i+1}"
        for i, coef in enumerate(coeffs)
    ])

def format_lp_problem(A, b, c, signs, sense):
    title = "МІНІМІЗАЦІЯ" if sense == 'min' else "МАКСИМІЗАЦІЯ"
    lines = [f"Пряма задача: {title}", f"F = {format_expression(c)}"]
    lines.append("Умови:")
    for row, sign, rhs in zip(A, signs, b):
        lines.append(f"  {format_expression(row)} {sign} {rhs}")
    return "\n".join(lines)

def format_dual_output(dual):
    title = "МІНІМІЗАЦІЯ" if dual['dual_sense'] == 'min' else "МАКСИМІЗАЦІЯ"
    lines = [f"Двоїста задача: {title}", f"Z = {format_expression(dual['c_dual'], vars_prefix='y')}"]
    lines.append("Умови:")
    for row, sign, rhs in zip(dual['A_dual'], dual['constraint_signs'], dual['b_dual']):
        lines.append(f"  {format_expression(row, vars_prefix='y')} {sign} {rhs}")
    return "\n".join(lines)

def init_table(c, A, b, signs, goal='min'):
    A = Matrix(A)
    b = Matrix(b).reshape(len(b), 1)
    vars = [f"y{i+1}" for i in range(len(c))]
    tab = A
    base = []
    arts = []

    for i, s in enumerate(signs):
        if s in ("<=", ">="):
            slack = Matrix.zeros(len(signs), 1)
            slack[i] = 1 if s == "<=" else -1
            tab = tab.row_join(slack)
            vars.append(f"y{len(vars)+1}")
            base.append(vars[-1])
            if s == ">=":
                art = Matrix.zeros(len(signs), 1)
                art[i] = 1
                tab = tab.row_join(art)
                art_name = f"w{len([v for v in vars if v.startswith('w')]) + 1}"
                vars.append(art_name)
                base[-1] = art_name
                arts.append((i, tab.cols - 1))

    tab = tab.row_join(b)

    cost = Matrix([0] * (tab.cols - 1) + [0])
    for i in range(len(c)):
        cost[i] = sp.Rational(c[i])

    sign = 1 if goal == 'min' else -1
    for idx, _ in arts:
        row = tab.row(idx)
        rhs = row[-1]
        for j in range(len(cost) - 1):
            cost[j] += sign * M * (-row[j])
        cost[-1] += sign * M * rhs

    for _, col in arts:
        cost[col] = 0

    for i in range(len(cost) - 1):
        cost[i] *= -1

    return tab, cost, vars, base

def solve_simplex(table, c, var_names, basis, goal='min'):
    num_vars = c.shape[0]
    num_constraints = table.shape[0]

    while True:
        numeric_c = [float(expr.subs(M, 1000)) for expr in c[:-1]]
        if goal == 'min':
            pivot_col = max((j for j in range(len(numeric_c)) if numeric_c[j] > 0), default=None, key=lambda j: numeric_c[j])
        else:
            pivot_col = min((j for j in range(len(numeric_c)) if numeric_c[j] < 0), default=None, key=lambda j: numeric_c[j])

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

        pivot_row = ratios.index(min(ratios))
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
        if basis[i] in var_names:
            solution[var_names.index(basis[i])] = table[i, -1]

    return solution, c[-1], basis

A = [ [1, 1], [3, 7], [1, 2], [1, 0], [0, 1]]
b = [8, 21, 6, 7.5, 7.5]
c = [1, 1]
signs = ['<=', '>=', '>=', '<=', '<=']
sense = 'min'

A_norm, b_norm, signs_norm = normalize_primal(A, b, signs, sense)
print(format_lp_problem(A_norm, b_norm, c, signs_norm, sense))
print()
dual = build_dual(A_norm, b_norm, c, sense)
dual['A_dual'], dual['b_dual'], dual['constraint_signs'] = fix_negative_rhs_for_dual(
    dual['A_dual'], dual['b_dual'], dual['constraint_signs']
)
print(format_dual_output(dual))
print("\n========== ОПТИМУМ ==========")
table, z_row, var_names, basis = init_table(
    dual['c_dual'], dual['A_dual'], dual['b_dual'],
    dual['constraint_signs'], goal=dual['dual_sense']
)
solution, optimal_value, final_basis = solve_simplex(
    table.copy(), z_row.copy(), var_names, basis.copy(), goal=dual['dual_sense']
)

if solution is not None:
    variables = [(var_names[i], simplify(solution[i])) for i in range(len(solution)) if simplify(solution[i]) != 0]
    print(tabulate(variables, headers=["Змінна", "Значення"], tablefmt="plain"))
    print(f"Оптимальне значення цільової функції: Z = {simplify(optimal_value)}")
else:
    print("Розв’язок не знайдено.")
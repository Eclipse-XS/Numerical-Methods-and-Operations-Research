import numpy as np
import sympy as sp
from sympy import Matrix, simplify

M = sp.Symbol("M")

def format_expr(coeffs, var='x'):
    return " ".join([f"{'' if a < 0 or i == 0 else '+'}{a}{var}{i+1}" for i, a in enumerate(coeffs)])

def fix_constraints(A, b, signs, goal, flip=False):
    target = '>=' if goal == 'min' else '<='
    new_A, new_b, new_signs = [], [], []
    for i in range(len(b)):
        row, val, sign = A[i], b[i], signs[i]
        if (flip and val < 0) or (not flip and sign != target):
            row, val = [-x for x in row], -val
            sign = {'>=': '<=', '<=': '>=', '=': '='}[sign]
        new_A.append(row)
        new_b.append(val)
        new_signs.append(sign)
    return new_A, new_b, new_signs

def build_dual(A, b, c, goal='min'):
    A, b, c = np.array(A), np.array(b), np.array(c)
    return {
        'type': 'max' if goal == 'min' else 'min',
        'A': A.T.tolist(),
        'b': c.tolist(),
        'c': b.tolist(),
        'signs': ['<=' if goal == 'min' else '>='] * len(A[0])
    }

def create_table(c, A, b, signs, goal='max'):
    n = len(c)
    A, b = Matrix(A), Matrix(b).reshape(len(b), 1)
    var_names = [f"x{i+1}" for i in range(n)]
    basis_vars, artificial_vars = [], []
    table = A
    next_var = n + 1

    for idx, cond in enumerate(signs):
        if cond == "<=" :
            col = Matrix.zeros(len(signs), 1)
            col[idx] = 1
            table = table.row_join(col)
            var_names.append(f"x{next_var}")
            basis_vars.append(f"x{next_var}")
            next_var += 1
        elif cond == ">=":
            neg_col = Matrix.zeros(len(signs), 1)
            neg_col[idx] = -1
            table = table.row_join(neg_col)
            var_names.append(f"x{next_var}")
            next_var += 1
            art_col = Matrix.zeros(len(signs), 1)
            art_col[idx] = 1
            table = table.row_join(art_col)
            art_name = f"w{len([v for v in var_names if v.startswith('w')])+1}"
            var_names.append(art_name)
            basis_vars.append(art_name)
            artificial_vars.append((idx, table.cols-1))
        elif cond == "=":
            art_col = Matrix.zeros(len(signs), 1)
            art_col[idx] = 1
            table = table.row_join(art_col)
            art_name = f"w{len([v for v in var_names if v.startswith('w')])+1}"
            var_names.append(art_name)
            basis_vars.append(art_name)
            artificial_vars.append((idx, table.cols-1))

    table = table.row_join(b)

    z = Matrix([0]*(table.cols-1) + [0])
    for i in range(len(c)):
        z[i] = sp.Rational(c[i])

    mult = 1 if goal == 'min' else -1
    for idx, _ in artificial_vars:
        for j in range(len(z)-1):
            z[j] += mult * M * (-table[idx, j])
        z[-1] += mult * M * table[idx, -1]

    for _, col in artificial_vars:
        z[col] = 0

    if goal == 'min':
        for j in range(len(z)-1):
            z[j] *= -1

    return table, z, var_names, basis_vars

def pick_pivot_col(z, goal='max'):
    z_vals = [float(expr.subs(M, 1000)) for expr in z[:-1]]
    if goal == 'max':
        return z_vals.index(min(z_vals)) if min(z_vals) < 0 else None
    else:
        pos = [i for i, v in enumerate(z_vals) if v > 0]
        return max(pos, key=lambda i: z_vals[i]) if pos else None

def simplex_method(table, z, var_names, basis_vars, goal='max'):
    while True:
        pivot_col = pick_pivot_col(z, goal)
        if pivot_col is None:
            break

        ratios = [(table[i, -1] / table[i, pivot_col] if table[i, pivot_col] > 0 else float('inf')) for i in range(table.rows)]
        if all(r == float('inf') for r in ratios):
            return None, None, None

        pivot_row = ratios.index(min(ratios))
        pivot_elem = table[pivot_row, pivot_col]

        table.row_op(pivot_row, lambda x, _: simplify(x / pivot_elem))
        for i in range(table.rows):
            if i != pivot_row:
                factor = table[i, pivot_col]
                table.row_op(i, lambda x, j: simplify(x - factor * table[pivot_row, j]))

        factor = z[pivot_col]
        z = Matrix([simplify(z[j] - factor * table[pivot_row, j]) for j in range(z.shape[0])])
        basis_vars[pivot_row] = var_names[pivot_col]

    sol = [0]*len(var_names)
    for i in range(table.rows):
        if basis_vars[i] in var_names[:len(c)]:
            sol[var_names.index(basis_vars[i])] = table[i, -1]

    return sol, z[-1], basis_vars

if __name__ == "__main__":
    A = [[1, 2], [2, 4], [6, 1], [1, 3]]
    b = [5, 16, 6, 9]
    c = [3, 1]
    signs = ['>=', '<=', '>=', '<=']
    goal = input("Enter the type of task (min/max): ").strip().lower()

    A_adj, b_adj, signs_adj = fix_constraints(A, b, signs, goal)
    print(f"Primal ({goal}):")
    print(f"Objective Function: F(x) = {format_expr(c)} -> {goal}")
    print("Constraints:")
    for i in range(len(b_adj)):
        print(f"  {format_expr(A_adj[i])} >= {b_adj[i]}")

    dual = build_dual(A_adj, b_adj, c, goal)
    print(f"\nDual ({dual['type']}):")
    print(f"Objective Function: Z(y) = {format_expr(dual['c'], 'y')} -> {dual['type']}")
    print("Constraints:")
    for i in range(len(dual['A'])):
        print(f"  {format_expr(dual['A'][i], 'y')} {dual['signs'][i]} {dual['b'][i]}")

    dual['A'], dual['b'], dual['signs'] = fix_constraints(dual['A'], dual['b'], dual['signs'], dual['type'], flip=True)

    table, z, var_names, basis_vars = create_table(dual['c'], dual['A'], dual['b'], dual['signs'], dual['type'])
    solution, optimum, final_basis = simplex_method(table, z, var_names, basis_vars, dual['type'])

    print("\nResult:")
    if solution:
        for i, val in enumerate(solution):
            if val != 0:
                print(f"y{i+1} = {simplify(val)}")
        print(f"Z(y) = {simplify(optimum)}")
    else:
        print("No feasible solution or problem is unbounded.")
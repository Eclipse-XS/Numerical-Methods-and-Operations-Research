import sympy as sp
from sympy import Matrix, simplify

M = sp.Symbol("M")

def render_ascii_table(table, names, bases, z_row, title=None):
    def row_text(cells, widths):
        return "| " + " | ".join(f"{str(c).center(w)}" for c, w in zip(cells, widths)) + " |"

    headers = ["Базис"] + [v for v in names] + ["B"]
    data = []

    for i in range(table.rows):
        row = [bases[i]] + [simplify(table[i, j]) for j in range(table.cols)]
        data.append(row)

    z_line = ["Z"] + [simplify(val) for val in z_row]
    data.append(z_line)

    widths = [max(len(str(item)) for item in col) for col in zip(*([headers] + data))]

    border = "+" + "+".join("-" * (w + 2) for w in widths) + "+"

    if title:
        print(f"\n{title.center(sum(widths) + 3 * len(widths))}\n")

    print(border)
    print(row_text(headers, widths))
    print(border)
    for row in data:
        print(row_text(row, widths))
        print(border)

def setup_initial_table(c, A, b, conditions, target='max'):
    n = len(c)
    A, b = Matrix(A), Matrix(b).reshape(len(b), 1)
    var_list = [f"x{i+1}" for i in range(n)]
    basis_list, artificial = [], []
    T = A
    next_idx = n + 1

    for idx, cond in enumerate(conditions):
        if cond == "<=":
            col = Matrix.zeros(len(conditions), 1)
            col[idx] = 1
            T = T.row_join(col)
            slack = f"x{next_idx}"
            var_list.append(slack)
            basis_list.append(slack)
            next_idx += 1
        elif cond == ">=":
            neg_col = Matrix.zeros(len(conditions), 1)
            neg_col[idx] = -1
            T = T.row_join(neg_col)
            surplus = f"x{next_idx}"
            var_list.append(surplus)
            next_idx += 1
            art_col = Matrix.zeros(len(conditions), 1)
            art_col[idx] = 1
            T = T.row_join(art_col)
            art_var = f"w{len([v for v in var_list if v.startswith('w')]) + 1}"
            var_list.append(art_var)
            basis_list.append(art_var)
            artificial.append((idx, T.cols - 1))
        elif cond == "=":
            art_col = Matrix.zeros(len(conditions), 1)
            art_col[idx] = 1
            T = T.row_join(art_col)
            art_var = f"w{len([v for v in var_list if v.startswith('w')]) + 1}"
            var_list.append(art_var)
            basis_list.append(art_var)
            artificial.append((idx, T.cols - 1))

    T = T.row_join(b)

    z = Matrix([0] * (T.cols - 1) + [0])
    for i in range(len(c)):
        z[i] = sp.Rational(c[i])

    mult = 1 if target == 'min' else -1
    for idx, _ in artificial:
        for j in range(len(z)-1):
            z[j] += mult * M * (-T[idx, j])
        z[-1] += mult * M * T[idx, -1]

    for _, col in artificial:
        z[col] = 0

    for j in range(len(z)-1):
        z[j] *= -1

    return T, z, var_list, basis_list

def find_pivot_col(z, goal='max'):
    z_num = [float(expr.subs(M, 1000)) for expr in z[:-1]]
    if goal == 'max':
        return z_num.index(min(z_num)) if min(z_num) < 0 else None
    if goal == 'min':
        pos = [i for i, v in enumerate(z_num) if v > 0]
        return max(pos, key=lambda i: z_num[i]) if pos else None

def run_simplex(table, z, var_list, basis_list, target='max'):
    render_ascii_table(table, var_list, basis_list, z, title="ПОЧАТКОВА ТАБЛИЦЯ")

    while True:
        pivot_col = find_pivot_col(z, target)
        if pivot_col is None:
            break

        ratios = [(table[i, -1] / table[i, pivot_col] if table[i, pivot_col] > 0 else float('inf')) for i in range(table.rows)]
        if all(r == float('inf') for r in ratios):
            return None, None, None

        pivot_row = ratios.index(min(ratios, key=float))
        pivot_elem = table[pivot_row, pivot_col]

        table.row_op(pivot_row, lambda x, _: simplify(x / pivot_elem))
        for i in range(table.rows):
            if i != pivot_row:
                factor = table[i, pivot_col]
                table.row_op(i, lambda x, j: simplify(x - factor * table[pivot_row, j]))

        factor = z[pivot_col]
        z = Matrix([simplify(z[j] - factor * table[pivot_row, j]) for j in range(z.shape[0])])
        basis_list[pivot_row] = var_list[pivot_col]

    render_ascii_table(table, var_list, basis_list, z, title="КІНЦЕВА ТАБЛИЦЯ")

    sol = [0] * len(var_list)
    for i in range(table.rows):
        if basis_list[i] in var_list[:len(c)]:
            sol[var_list.index(basis_list[i])] = table[i, -1]

    return sol, z[-1], basis_list

# === INPUT DATA ===
A = [[1, 2], [2, 4], [6, 1], [1, 3]]
b = [5, 16, 6, 9]
c = [3, 1]
signs = ['>=', '<=', '>=', '<=']
goal = 'min'

table, z_row, var_names, basis_vars = setup_initial_table(c, A, b, signs, goal)
solution, optimal_value, final_basis = run_simplex(table.copy(), z_row.copy(), var_names, basis_vars.copy(), goal)

# === OUTPUT ===
x_vals = [simplify(solution[i]) for i in range(len(c))]
print("\nОптимальне рішення задачі:")
print(f"(x₁ = {x_vals[0]}; x₂ = {x_vals[1]})")
print(f"Значення цільової функції: F(x) = {simplify(optimal_value)}")
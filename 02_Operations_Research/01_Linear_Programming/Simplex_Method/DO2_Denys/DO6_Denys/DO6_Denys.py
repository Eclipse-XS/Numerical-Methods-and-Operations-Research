import sympy as sp
from sympy import Matrix, simplify
from tabulate import tabulate

M = sp.Symbol("M")

# Кольори
YELLOW = "\033[43m"
GREEN = "\033[42m"
BLUE = "\033[44m"
RESET = "\033[0m"

def build_canonical_table(c, A, b, signs, goal='max'):
    num_vars = len(c)
    A = Matrix(A)
    b = Matrix(b).reshape(len(b), 1)
    var_names = [f"x{i+1}" for i in range(num_vars)]
    basis = []
    artificial_vars = []

    table = A

    next_var_idx = num_vars + 1

    for i, sign in enumerate(signs):
        if sign == "<=":
            slack_col = Matrix.zeros(len(signs), 1)
            slack_col[i] = 1
            table = table.row_join(slack_col)
            slack_var = f"x{next_var_idx}"
            next_var_idx += 1
            var_names.append(slack_var)
            basis.append(slack_var)

        elif sign == ">=":
            # surplus
            surplus_col = Matrix.zeros(len(signs), 1)
            surplus_col[i] = -1
            table = table.row_join(surplus_col)
            surplus_var = f"x{next_var_idx}"
            next_var_idx += 1
            var_names.append(surplus_var)

            # artificial
            artificial_col = Matrix.zeros(len(signs), 1)
            artificial_col[i] = 1
            table = table.row_join(artificial_col)
            art_var = f"a{len([v for v in var_names if v.startswith('w')]) + 1}"
            var_names.append(art_var)
            basis.append(art_var)
            artificial_vars.append((i, table.cols - 1))

        elif sign == "=":
            # artificial only
            artificial_col = Matrix.zeros(len(signs), 1)
            artificial_col[i] = 1
            table = table.row_join(artificial_col)
            art_var = f"w{len([v for v in var_names if v.startswith('w')]) + 1}"
            var_names.append(art_var)
            basis.append(art_var)
            artificial_vars.append((i, table.cols - 1))

    # Add RHS
    table = table.row_join(b)

    # Z row (cost row)
    cost_row = Matrix([0] * (table.cols - 1) + [0])
    for i in range(len(c)):
        cost_row[i] = sp.Rational(c[i])

    sign_mult = 1 if goal == 'min' else -1
    for row_idx, _ in artificial_vars:
        row = table.row(row_idx)
        rhs = row[-1]
        for j in range(len(cost_row) - 1):
            cost_row[j] += sign_mult * M * (-row[j])
        cost_row[-1] += sign_mult * M * rhs

    for _, art_col in artificial_vars:
        cost_row[art_col] = 0

    for i in range(len(cost_row) - 1):
        cost_row[i] *= -1

    return table, cost_row, var_names, basis



def print_colored_table(table, var_names, basis, z_row, pivot_row=None, pivot_col=None, title=None):
    headers = ["Базис"] + [v.upper() for v in var_names] + ["B"]
    data = []

    for i in range(table.rows):
        row_data = []
        for j in range(table.cols):
            val = str(simplify(table[i, j]))
            if i == pivot_row and j == pivot_col:
                val = BLUE + val + RESET
            elif i == pivot_row:
                val = GREEN + val + RESET
            elif j == pivot_col:
                val = YELLOW + val + RESET
            row_data.append(val)
        data.append([basis[i]] + row_data)

    z_data = []
    for j in range(len(z_row)):
        val = str(simplify(z_row[j]))
        if j == pivot_col:
            val = YELLOW + val + RESET
        z_data.append(val)

    data.append(["F"] + z_data)

    if title:
        print(f"\n📌 {title}")
    print(tabulate(data, headers=headers, tablefmt="grid", stralign="center"))


def choose_pivot_column(c, goal='max'):
    numeric_c = [float(expr.subs(M, 1000)) for expr in c[:-1]]
    if goal == 'max':
        return numeric_c.index(min(numeric_c)) if min(numeric_c) < 0 else None
    elif goal == 'min':
        pos = [i for i, val in enumerate(numeric_c) if val > 0]
        return max(pos, key=lambda i: numeric_c[i]) if pos else None


def symbolic_simplex_verbose(table, c, var_names, basis, goal='max'):
    num_vars = c.shape[0]
    num_constraints = table.shape[0]

    iteration = 1

    print_colored_table(table, var_names, basis, c, title="ПОЧАТКОВА ТАБЛИЦЯ")

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
            print("❌ Розв’язок не обмежений.")
            return None, None, None

        pivot_row = ratios.index(min(ratios, key=lambda r: float(r)))
        pivot_element = table[pivot_row, pivot_col]

        print_colored_table(table, var_names, basis, c, pivot_row, pivot_col, title=f"ІТЕРАЦІЯ {iteration}")

        table.row_op(pivot_row, lambda x, _: simplify(x / pivot_element))

        for i in range(num_constraints):
            if i != pivot_row:
                factor = table[i, pivot_col]
                table.row_op(i, lambda x, j: simplify(x - factor * table[pivot_row, j]))

        factor = c[pivot_col]
        c = Matrix([simplify(c[j] - factor * table[pivot_row, j]) for j in range(num_vars)])
        basis[pivot_row] = var_names[pivot_col]
        iteration += 1

    print_colored_table(table, var_names, basis, c, title="КІНЦЕВА ТАБЛИЦЯ")

    solution = [0] * num_vars
    for i in range(num_constraints):
        col = var_names.index(basis[i])
        solution[col] = table[i, -1]

    return solution, c[-1], basis
A = [
    [8,  -5],
    [2,   5],
    [-6,  5],
    [2,   1]
]

b = [40, 10, 60, 14]

signs = ['<=', '>=', '<=', '<=']

c = [3, 4]  # F(x1, x2) = 3x1 + 4x2

goal = 'max'  # або 'max' — залежно від варіанту, вибирається один


table, z_row, var_names, basis_vars = build_canonical_table(c, A, b, signs, goal)
solution, optimal_value, final_basis = symbolic_simplex_verbose(table.copy(), z_row.copy(), var_names, basis_vars.copy(), goal)

# === РЕЗУЛЬТАТ ===
print("\n========== РЕЗУЛЬТАТ ==========")
for i, val in enumerate(solution):
    if val != 0:
        print(f"{var_names[i]} = {simplify(val)}")
print(f"Оптимальне значення цільової функції: {simplify(optimal_value)}")

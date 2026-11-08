import sympy as sp
from sympy import Matrix, simplify
from tabulate import tabulate

M = sp.Symbol("M")

YELLOW = "\033[43m"
GREEN = "\033[42m"
BLUE = "\033[44m"
RESET = "\033[0m"

def convert_fractional_to_linear_lp(c_num, c_den, A, b, signs, goal):
    c_num = list(map(sp.Rational, c_num))
    c_den = list(map(sp.Rational, c_den))
    A = [[sp.Rational(aij) for aij in row] for row in A]
    b = list(map(sp.Rational, b))

    print("Дробово-лінійна задача:")
    x_vars = [sp.Symbol(f"x{i+1}") for i in range(len(c_num))]
    F_expr = sp.Add(*[a * x for a, x in zip(c_num, x_vars)]) / sp.Add(*[a * x for a, x in zip(c_den, x_vars)])
    print(f"F = {F_expr} --> {goal}")
    print("За умов:")
    for i in range(len(A)):
        lhs = sp.Add(*[A[i][j] * x_vars[j] for j in range(len(x_vars))])
        print(lhs, signs[i], b[i])
    print()

    n = len(c_num)
    c_lin = list(c_num) + [0]

    A_lin, new_signs, b_lin = [], [], []
    for i in range(len(A)):
        row = list(A[i]) + [-b[i]]
        if signs[i] == ">=":
            row = [-a for a in row]
            new_signs.append("<=")
            b_lin.append(0)
        elif signs[i] == "<=":
            new_signs.append("<=")
            b_lin.append(0)
        elif signs[i] == "=":
            new_signs.append("=")
            b_lin.append(0)
        A_lin.append(row)

    denom_constraint = list(c_den) + [0]
    A_lin.append(denom_constraint)
    new_signs.append("=")
    b_lin.append(1)

    columns = [f"y{i+1}" for i in range(n)] + ["y0"]
    print("Лінійна задача:")
    y_vars = [sp.Symbol(name) for name in columns]
    L_expr = sp.Add(*[a * y for a, y in zip(c_lin, y_vars)])
    print(f"f = {L_expr} --> {goal}")
    print("За умов:")
    for i in range(len(A_lin)):
        lhs = sp.Add(*[A_lin[i][j] * y_vars[j] for j in range(len(columns))])
        print(lhs, new_signs[i], b_lin[i])
    print()

    return A_lin, c_lin, columns, new_signs, b_lin

def build_canonical_table(c, A, b, signs, goal='max'):
    A = Matrix(A)
    b = Matrix(b).reshape(len(b), 1)
    table = A
    var_names = [f"y{i+1}" for i in range(A.shape[1] - 1)] + ["y0"]
    basis = []
    artificial_vars = []

    slack_surplus_count = 0
    artificial_count = 1

    for i, sign in enumerate(signs):
        if sign == "<=":
            slack = f"s{slack_surplus_count+1}"
            slack_col = Matrix.eye(len(signs)).col(i)
            table = table.row_join(slack_col)
            var_names.append(slack)
            basis.append(slack)
            slack_surplus_count += 1
        elif sign == ">=":
            surplus = f"s{slack_surplus_count+1}"
            surplus_col = -Matrix.eye(len(signs)).col(i)
            artificial = f"w{artificial_count}"
            artificial_col = Matrix.eye(len(signs)).col(i)
            table = table.row_join(surplus_col)
            table = table.row_join(artificial_col)
            var_names.extend([surplus, artificial])
            basis.append(artificial)
            artificial_vars.append((i, artificial))
            slack_surplus_count += 1
            artificial_count += 1
        elif sign == "=":
            artificial = f"a{artificial_count}"
            artificial_col = Matrix.eye(len(signs)).col(i)
            table = table.row_join(artificial_col)
            var_names.append(artificial)
            basis.append(artificial)
            artificial_vars.append((i, artificial))
            artificial_count += 1

    table = table.row_join(b)

    cost_row = [0] * (table.cols - 1)
    for i in range(len(c)):
        cost_row[i] = sp.Rational(c[i])
    cost_row = Matrix(cost_row + [0])

    sign_mult = 1 if goal == 'min' else -1
    for idx, a_var in artificial_vars:
        rhs = table[idx, -1]
        for j in range(table.cols - 1):
            if var_names[j] != a_var:
                cost_row[j] += sign_mult * M * (-table[idx, j])
        cost_row[var_names.index(a_var)] = 0
        cost_row[-1] += sign_mult * M * rhs

    for i in range(len(cost_row) - 1):
        cost_row[i] *= -1

    return table, cost_row, var_names, basis

def print_colored_table(table, var_names, basis, z_row, pivot_row=None, pivot_col=None, title=None):
    headers = ["Базис"] + [v.upper() for v in var_names] + ["B"]
    data = []

    def format_cell(cell, row_idx, col_idx):
        cell = simplify(cell)
        val = str(cell)
        if row_idx == pivot_row and col_idx == pivot_col:
            return f"{BLUE}{val}{RESET}"
        elif row_idx == pivot_row:
            return f"{GREEN}{val}{RESET}"
        elif col_idx == pivot_col:
            return f"{YELLOW}{val}{RESET}"
        else:
            return val

    for i in range(table.rows):
        row = [basis[i]] + [format_cell(table[i, j], i, j) for j in range(table.cols)]
        data.append(row)

    z = ["Z"] + [format_cell(z_row[j], -1, j) for j in range(len(z_row))]
    data.append(z)

    if title:
        print(f"\n📌 {title}")
    print(tabulate(data, headers=headers, tablefmt="fancy_grid", stralign="center"))


def symbolic_simplex_verbose(table, c, var_names, basis, goal='max'):
    num_vars = c.shape[0]
    num_constraints = table.shape[0]
    iteration = 0

    print_colored_table(table, var_names, basis, c, title="ПОЧАТКОВА СИМПЛЕКС-ТАБЛИЦЯ")

    while True:
        numeric_c = [float(expr.subs(M, 1000)) for expr in c[:-1]]
        if goal == 'max':
            pivot_col = numeric_c.index(min(numeric_c)) if min(numeric_c) < 0 else None
        else:
            pos = [i for i, val in enumerate(numeric_c) if val > 0]
            pivot_col = max(pos, key=lambda i: numeric_c[i]) if pos else None

        if pivot_col is None:
            print_colored_table(table, var_names, basis, c, title="КІНЦЕВА СИМПЛЕКС-ТАБЛИЦЯ")
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

        # 🔹 Покажемо таблицю ДО обчислення (але вже знаємо, що оберемо)
        print_colored_table(table, var_names, basis, c, pivot_row, pivot_col, title=f"ІТЕРАЦІЯ {iteration}")

        # 🔧 Обчислення
        pivot_element = table[pivot_row, pivot_col]
        table.row_op(pivot_row, lambda x, _: simplify(x / pivot_element))
        for i in range(num_constraints):
            if i != pivot_row:
                factor = table[i, pivot_col]
                table.row_op(i, lambda x, j: simplify(x - factor * table[pivot_row, j]))

        factor = c[pivot_col]
        c = Matrix([simplify(c[j] - factor * table[pivot_row, j]) for j in range(num_vars)])
        basis[pivot_row] = var_names[pivot_col]
        iteration += 1

    solution = [0] * num_vars
    for i in range(num_constraints):
        col = var_names.index(basis[i])
        solution[col] = table[i, -1]

    return solution, c[-1], basis

# === Вхідні дані ===
goal = 'max'  # або 'max', якщо потрібно максимізувати

# Чисельник і знаменник цільової функції:
c_num = [2, 2]      # чисельник: 2x1 + 2x2
c_den = [2, 3]      # знаменник: 2x1 + 3x2

# Матриця коефіцієнтів обмежень:
A = [
    [-3, 2],      # 3x1 - 2x2 ≥ -6
    [1,  2],      # x1 + 2x2 ≥ 3
    [1,  0],      # x1 ≤ 3
    [0,  1],      # x2 ≤ 5
]

# Вільні члени:
b = [6, 3, 3, 5]

# Знаки обмежень:
signs = ["<=", ">=", "<=", "<="]




# === Виконання ===
A_lin, c_lin, columns, signs_new, b_lin = convert_fractional_to_linear_lp(c_num, c_den, A, b, signs, goal)
table, cost_row, var_names, basis = build_canonical_table(c_lin, A_lin, b_lin, signs_new, goal)

# === Симплекс-метод ===
solution, opt_value, final_basis = symbolic_simplex_verbose(table.copy(), cost_row.copy(), var_names, basis, goal)

# === Відображення результатів ===
print(f"\nОптимум лінійної функції G(y): {opt_value}")

print("Зворотня заміна:")
y0_index = var_names.index("y0")
y0_val = solution[y0_index]

if y0_val == 0:
    print("y0 = 0 → неможливо виконати зворотну заміну (ділення на 0)")
else:
    x_vals = []
    for i in range(len(c_num)):
        yi = solution[var_names.index(f"y{i+1}")]
        xi = simplify(yi / y0_val)
        x_vals.append(xi)
        print(f"x{i+1} = y{i+1}/y0 = {yi} / {y0_val} = {xi}")

    numerator = sum(c_num[i] * x_vals[i] for i in range(len(c_num)))
    denominator = sum(c_den[i] * x_vals[i] for i in range(len(c_den)))
    F_val = simplify(numerator / denominator)
    print(f"Оптимум дробово-раціональної функції F(x): {F_val}")

import numpy as np
import pandas as pd
import sympy as sp

def convert_fractional_to_linear_lp(c_num, c_den, A, b):
    c_num = np.array(c_num, dtype=float)
    c_den = np.array(c_den, dtype=float)
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)

    print("Дробово-лінійна задача:")
    x_vars = [sp.Symbol(f"x{i+1}") for i in range(len(c_num))]
    F_expr = sp.Add(*[a * x for a, x in zip(c_num, x_vars)]) / sp.Add(*[a * x for a, x in zip(c_den, x_vars)])
    print("F(x) =", F_expr)
    print("Система обмежень:")
    for i in range(len(A)):
        lhs = sp.Add(*[A[i][j] * x_vars[j] for j in range(len(x_vars))])
        sign = "<=" if b[i] >= 0 else ">="
        rhs = abs(b[i])
        print(lhs, sign, rhs)
    print()

    n = len(c_num)
    m = len(A)
    c_lin = list(c_num) + [0]

    A_lin, signs, b_lin = [], [], []
    for i in range(m):
        if b[i] >= 0:
            row = list(A[i]) + [-b[i]]
            signs.append("<=")
        else:
            row = [-a for a in A[i]] + [abs(b[i])]
            signs.append("<=")
        A_lin.append(row)
        b_lin.append(0)

    denom_constraint = list(c_den) + [0]
    A_lin.append(denom_constraint)
    signs.append("=")
    b_lin.append(1)

    columns = [f"y{i+1}" for i in range(n)] + ["y0"]
    df = pd.DataFrame(A_lin, columns=columns)
    df["sign"] = signs
    df["b"] = b_lin

    print("Лінійна задача:")
    y_vars = [sp.Symbol(name) for name in columns]
    L_expr = sp.Add(*[a * y for a, y in zip(c_lin, y_vars)])
    print("L(y) =", L_expr)
    print("Система обмежень:")
    for i in range(len(df)):
        lhs = sp.Add(*[df.iloc[i][col] * sp.Symbol(col) for col in columns])
        sign = df.iloc[i]['sign']
        rhs = df.iloc[i]['b']
        print(lhs, sign, rhs)
    print()

    return df, c_lin, x_vars, c_num

def to_canonical_form(df, c_lin, goal='min'):
    canonical_df = df.copy()
    slack_vars, artificial_vars = [], []
    M = sp.Symbol("M")
    artificial_count = 1
    original_vars = [col for col in df.columns if col.startswith("y")]

    for i, sign in enumerate(df["sign"]):
        if sign == "<=" :
            slack = f"y{len(original_vars) + len(slack_vars)}"
            canonical_df[slack] = 0
            canonical_df.loc[i, slack] = 1
            slack_vars.append(slack)
        elif sign == ">=":
            surplus = f"y{len(original_vars) + len(slack_vars) + 1}"
            artificial = f"w{artificial_count}"
            artificial_count += 1
            canonical_df[surplus] = 0
            canonical_df[artificial] = 0
            canonical_df.loc[i, surplus] = -1
            canonical_df.loc[i, artificial] = 1
            slack_vars.append(surplus)
            artificial_vars.append(artificial)
        elif sign == "=":
            artificial = f"w{artificial_count}"
            artificial_count += 1
            canonical_df[artificial] = 0
            canonical_df.loc[i, artificial] = 1
            artificial_vars.append(artificial)

    canonical_df.drop(columns=["sign"], inplace=True)
    vars_in_table = [col for col in canonical_df.columns if col != "b"]

    c_sym = []
    for col in vars_in_table:
        if col.startswith("y") and col != "y0":
            index = int(col[1:]) - 1
            if index < len(c_lin) - 1:
                c_sym.append(sp.Rational(c_lin[index]))
            else:
                c_sym.append(0)
        elif col == "y0":
            c_sym.append(sp.Rational(c_lin[-1]))
        else:
            c_sym.append(0)
    c_sym.append(0)

    sign = 1 if goal == 'min' else -1
    for a_var in artificial_vars:
        idx = canonical_df.index[canonical_df[a_var] == 1][0]
        rhs = canonical_df.loc[idx, 'b']
        for i, col in enumerate(vars_in_table):
            if col != a_var:
                c_sym[i] += sign * M * (-canonical_df.loc[idx, col])
        c_sym[vars_in_table.index(a_var)] = 0
        c_sym[-1] += sign * M * rhs

    return canonical_df, c_sym, vars_in_table, M

def print_constraints_matrix_only(df):
    print("Канонічна форма:")
    vars_in_table = [col for col in df.columns if col != "b"]
    matrix = df[vars_in_table + ['b']]
    print(matrix.to_string(index=False))

def choose_pivot_column(c, goal):
    if goal == 'max':
        return np.argmin(c[:-1]) if np.min(c[:-1]) < 0 else None
    elif goal == 'min':
        pos = [j for j in range(len(c) - 1) if c[j] > 0]
        return max(pos, key=lambda j: c[j]) if pos else None

def simplex_method(table, c, goal):
    num_vars = len(c) - 1
    num_constraints = len(table)
    basis = []

    for col in range(num_vars):
        column = table[:, col]
        if list(column).count(1) == 1 and list(column).count(0) == num_constraints - 1:
            row = list(column).index(1)
            if len(basis) < num_constraints:
                basis.append(col)

    while True:
        pivot_col = choose_pivot_column(c, goal)
        if pivot_col is None:
            break

        ratios = [
            table[i, -1] / table[i, pivot_col] if table[i, pivot_col] > 0 else np.inf
            for i in range(num_constraints)
        ]
        if all(r == np.inf for r in ratios):
            print("Задача не має обмеженого розв’язку.")
            return None, None, None

        pivot_row = np.argmin(ratios)
        basis[pivot_row] = pivot_col
        table[pivot_row, :] /= table[pivot_row, pivot_col]

        for i in range(num_constraints):
            if i != pivot_row:
                table[i, :] -= table[i, pivot_col] * table[pivot_row, :]

        c -= c[pivot_col] * table[pivot_row, :]

    solution = np.zeros(num_vars)
    for i in range(num_constraints):
        solution[basis[i]] = table[i, -1]

    return solution, c[-1], basis



c_num = [2, 2]          # Чисельник цільової функції: 2x1 + 2x2
c_den = [-2, 3]         # Знаменник цільової функції: -2x1 + 3x2

A = [
    [5, -2],            # 5x1 - 2x2 ≤ 7
    [-1, 1],            # -x1 + x2 ≤ 5
    [1, 1],             # x1 + x2 ≤ 6
]

b = [7, 5, 6]

df, c_lin, x_vars, c_num = convert_fractional_to_linear_lp(c_num, c_den, A, b)
canonical_df_min, c_sym_min, vars_in_table, M = to_canonical_form(df, c_lin, goal='min')
canonical_df_max, c_sym_max, _, _ = to_canonical_form(df, c_lin, goal='max')

print_constraints_matrix_only(canonical_df_min)

z_min = [sp.simplify(c) for c in c_sym_min]
z_max = [sp.simplify(c) for c in c_sym_max]
z_row_min = [-c for c in z_min[:-1]] + [z_min[-1]]
z_row_max = [-c for c in z_max[:-1]] + [z_max[-1]]
print("\nОціночний рядок (min):", z_row_min)
print("Оціночний рядок (max):", z_row_max)

M_value = 1e6

print("\n=======МІНІМІЗАЦІЯ=======")
table_min = canonical_df_min.to_numpy(dtype=float)
z_row_min_np = np.array([float(expr.subs(M, M_value).evalf()) for expr in z_row_min])
sol_min, val_min, _ = simplex_method(table_min, z_row_min_np.copy(), goal='min')
sol_dict_min = {v: x for v, x in zip(vars_in_table, sol_min)}
y0_min = sol_dict_min.get("y0", 1)
x_result_min = [(f"x{i+1}", sol_dict_min.get(f"y{i+1}", 0) / y0_min if y0_min != 0 else float('inf')) for i in range(len(c_num))]
print("Оптимальний розв'язок:")
print("( " + ", ".join([f"y{i+1}={sol_dict_min.get(f'y{i+1}', 0):.3f}" for i in range(len(c_num))]) + " )")
print("L(y) =", round(val_min, 3))
print("Обернена заміна: yi -> xi :")
print("( " + ", ".join([f"{var}={val:.3f}" for var, val in x_result_min]) + " )")
c_num_expr = sum(c_num[i] * x_vars[i] for i in range(len(c_num)))
print("F(x) =", round(float(c_num_expr.subs({x_vars[i]: val for i, (_, val) in enumerate(x_result_min)}).evalf()), 2))

print("\n=======МАКСИМІЗАЦІЯ=======")
table_max = canonical_df_max.to_numpy(dtype=float)
z_row_max_np = np.array([float(expr.subs(M, M_value).evalf()) for expr in z_row_max])
sol_max, val_max, _ = simplex_method(table_max, z_row_max_np.copy(), goal='max')
sol_dict_max = {v: x for v, x in zip(vars_in_table, sol_max)}
y0_max = sol_dict_max.get("y0", 1)
x_result_max = [(f"x{i+1}", sol_dict_max.get(f"y{i+1}", 0) / y0_max if y0_max != 0 else float('inf')) for i in range(len(c_num))]
print("Оптимальний розв'язок:")
print("( " + ", ".join([f"y{i+1}={sol_dict_max.get(f'y{i+1}', 0):.3f}" for i in range(len(c_num))]) + " )")
print("L(y) =", round(val_max, 3))
print("Обернена заміна: yi -> xi :")
print("( " + ", ".join([f"{var}={val:.3f}" for var, val in x_result_max]) + " )")
print("F(x) =", float(c_num_expr.subs({x_vars[i]: val for i, (_, val) in enumerate(x_result_max)}).evalf()))

print("\n" + "=" * 40)
print("Група: ОІ-23")
print("Виконав: Гончаров Вадим Олегович")
print("Номер ЛР: №6")
print("Варіант: №75")
print("=" * 40)
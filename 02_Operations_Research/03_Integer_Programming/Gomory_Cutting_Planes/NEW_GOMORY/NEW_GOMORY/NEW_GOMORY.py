import sympy as sp
import numpy as np
from sympy import Matrix, nsimplify, simplify
from tabulate import tabulate
from termcolor import colored
M = sp.Symbol("M")

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
            surplus_col = Matrix.zeros(len(signs), 1)
            surplus_col[i] = -1
            table = table.row_join(surplus_col)
            surplus_var = f"x{next_var_idx}"
            next_var_idx += 1
            var_names.append(surplus_var)

            artificial_col = Matrix.zeros(len(signs), 1)
            artificial_col[i] = 1
            table = table.row_join(artificial_col)
            art_var = f"a{len([v for v in var_names if v.startswith('w')]) + 1}"
            var_names.append(art_var)
            basis.append(art_var)
            artificial_vars.append((i, table.cols - 1))

        elif sign == "=":
            artificial_col = Matrix.zeros(len(signs), 1)
            artificial_col[i] = 1
            table = table.row_join(artificial_col)
            art_var = f"a{len([v for v in var_names if v.startswith('w')]) + 1}"
            var_names.append(art_var)
            basis.append(art_var)
            artificial_vars.append((i, table.cols - 1))

    table = table.row_join(b)

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

def choose_pivot_column(c, goal='max'):
    numeric_c = [float(sp.sympify(expr).subs(M, 1000)) for expr in c[:-1]]
    if goal == 'max':
        return numeric_c.index(min(numeric_c)) if min(numeric_c) < 0 else None
    elif goal == 'min':
        pos = [i for i, val in enumerate(numeric_c) if val > 0]
        return max(pos, key=lambda i: numeric_c[i]) if pos else None

def symbolic_simplex_method_final(table, c, var_names, basis, goal='max'):
    table = np.array(table.tolist(), dtype=object)
    c = [sp.sympify(val) for val in c]
    num_constraints = table.shape[0]
    num_vars = table.shape[1] - 1

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
            return None, None, None, None

        pivot_row = ratios.index(min(ratios, key=lambda r: float(r)))
        pivot_element = table[pivot_row, pivot_col]

        table[pivot_row] = [simplify(x / pivot_element) for x in table[pivot_row]]

        for i in range(num_constraints):
            if i != pivot_row:
                factor = table[i, pivot_col]
                table[i] = [simplify(table[i][j] - factor * table[pivot_row][j]) for j in range(table.shape[1])]

        factor = c[pivot_col]
        c = [simplify(c[j] - factor * table[pivot_row][j]) for j in range(len(c))]
        basis[pivot_row] = var_names[pivot_col]

    solution = [0] * num_vars
    for i in range(num_constraints):
        if basis[i] in var_names[:num_vars]:
            idx = var_names.index(basis[i])
            solution[idx] = table[i, -1]

    return table, c, basis, solution

def is_integer(x): return np.isclose(float(x), round(float(x)))
def get_fractional_part(x): return float(x) - np.floor(float(x))

def find_gomory_row_index(table, var_names, basis):
    target_vars = ["x1", "x2"]
    row_candidates = []

    for i, base_var in enumerate(basis):
        if base_var in target_vars:
            frac_part = get_fractional_part(table[i, -1]) if not is_integer(table[i, -1]) else 0
            row_candidates.append((i, frac_part))

    if not row_candidates:
        return None

    max_row = max(row_candidates, key=lambda x: x[1])

    if max_row[1] == 0:
        return None

    return max_row[0]

def add_gomory_cut(table, c, var_names, basis):
    row_index = find_gomory_row_index(table, var_names, basis)
    if row_index is None:
        return table, c, var_names, basis

    plan_vars = []
    for var in ["x1", "x2"]:
        if var in basis:
            plan_vars.append(var)
        else:
            plan_vars.append("0")
    plan_str = "; ".join(plan_vars)
    print(f"\nОпорний план: ({plan_str})")

    rhs_val = float(table[row_index, -1])
    frac_value = Fraction(rhs_val - np.floor(rhs_val)).limit_denominator(1000)
    print(f"Максимальна дробова частина у {basis[row_index]}: β = {frac_value}")

    row = table[row_index]

    gomory_coeffs = [-(get_fractional_part(row[i])) for i in range(len(row) - 1)]
    gomory_rhs = -get_fractional_part(row[-1])

    print("\nРахуємо нові коефіцієнти α:")
    alpha_count = 1
    for idx, val in enumerate(row[:-1]):
        val_float = float(val)
        frac_val = Fraction(val_float).limit_denominator(1000)
        floor_val = int(np.floor(val_float))
        frac_part = frac_val - floor_val
        if np.isclose(frac_part, 0):
            continue
        print(f"α{alpha_count} = {{ {frac_val} }} = {frac_val} - ⌊{frac_val}⌋ = {frac_val} - ({floor_val}) = {frac_part}")
        alpha_count += 1

    cut_expr = ""
    for coeff, var in zip(gomory_coeffs, var_names):
        if not np.isclose(coeff, 0):
            coeff_frac = Fraction(coeff).limit_denominator(1000)
            term = f"{coeff_frac}*{var}" if coeff_frac >= 0 else f"- {abs(coeff_frac)}*{var}"
            cut_expr += (" + " if cut_expr and coeff_frac >= 0 else "") + term

    # 🔧 Заміна нумерації cut-змінних
    cut_count = sum(1 for v in var_names if v.startswith("cut"))
    new_var = f"cut{cut_count + 1}"

    cut_expr += f" + 1*{new_var}"

    rhs_frac = Fraction(gomory_rhs).limit_denominator(1000)

    print(f"\nДодано обмеження: {cut_expr.strip()} = {rhs_frac}")

    gomory_row = [Fraction(-get_fractional_part(row[i])).limit_denominator(1000) for i in range(len(row) - 1)]
    gomory_row += [1, Fraction(-get_fractional_part(row[-1])).limit_denominator(1000)]

    new_table = np.hstack((table[:, :-1], np.zeros((table.shape[0], 1)), table[:, -1:]))
    new_table = np.vstack((new_table, np.array(gomory_row, dtype=object)))
    new_c = np.append(c[:-1], [0, c[-1]])
    var_names.append(new_var)
    basis.append(new_var)

    return new_table, new_c, var_names, basis



from fractions import Fraction

def format_cell(val, eps=1e-10, max_denominator=1000):
    val = sp.nsimplify(val, rational=True)
    
    if val.free_symbols:
        return str(val)

    try:
        float_val = float(val)
        if abs(float_val) < eps:
            return "0"
        frac = Fraction(float_val).limit_denominator(max_denominator)
        return str(frac)
    except:
        return str(val)

def print_table_with_tabulate(table, var_names, basis, z_row):
    headers = ["Базис"] + [v.upper() for v in var_names] + ["B"]
    data = []
    for i in range(table.rows):
        row = [basis[i]] + [format_cell(cell) for cell in table.row(i)]
        data.append(row)
    z = ["Z"] + [format_cell(e) for e in z_row]
    data.append(z)
    print(tabulate(data, headers=headers, tablefmt="grid", stralign="center"))
def print_highlighted_table(table, c, var_names, basis, pivot_row, pivot_col):
    headers = ["Базис"] + [v.upper() for v in var_names] + ["B"]
    data = []

    for i, row in enumerate(table):
        row_data = []
        for j, val in enumerate(row):
            val_str = format_cell(val)
            if j == pivot_col and i == pivot_row:
                val_str = colored(val_str, "green", attrs=["reverse", "bold"])
            elif j == pivot_col:
                val_str = colored(val_str, "yellow")
            elif i == pivot_row:
                val_str = colored(val_str, "green")
            row_data.append(val_str)
        data.append([basis[i]] + row_data)

    z_row = []
    for j, val in enumerate(c):
        val_str = format_cell(val)
        if j == pivot_col:
            val_str = colored(val_str, "yellow")
        z_row.append(val_str)
    data.append(["Z"] + z_row)

    print(tabulate(data, headers=headers, tablefmt="grid", stralign="center"))


def dual_simplex(table, c, var_names, basis):
    num_constraints, num_vars = table.shape[0], table.shape[1] - 1
    iter_count = 1

    while True:
        if all(table[i, -1] >= 0 for i in range(num_constraints)):
            break

        pivot_row = min((i for i in range(num_constraints) if table[i, -1] < 0), key=lambda i: table[i, -1])

        row = table[pivot_row, :-1]
        candidates = []
        for j in range(num_vars):
            if row[j] < 0:
                ratio = abs(c[j] / row[j])
                candidates.append((ratio, j))

        if not candidates:
            print("❗ Немає допустимого опорного стовпця")
            return None, None, None, None

        _, pivot_col = min(candidates, key=lambda x: float(x[0].subs(M, 1000) if isinstance(x[0], sp.Basic) else x[0]))

        pivot_element = table[pivot_row, pivot_col]

        print(f"\n=== ІТЕРАЦІЯ {iter_count} ГОМОРІ ===")
        print(f"🔹 Опорний стовпець: {pivot_col + 1}, 🔹 Опорний рядок: {pivot_row + 1}, значення: {nsimplify(pivot_element, rational=True)}")
        print_highlighted_table(table, c, var_names, basis, pivot_row, pivot_col)
        iter_count += 1

        table[pivot_row] = [simplify(x / pivot_element) for x in table[pivot_row]]

        for i in range(num_constraints):
            if i != pivot_row:
                factor = table[i, pivot_col]
                table[i] = [simplify(table[i][j] - factor * table[pivot_row][j]) for j in range(table.shape[1])]

        factor = c[pivot_col]
        c = [simplify(c[j] - factor * table[pivot_row][j]) for j in range(len(c))]

        basis[pivot_row] = var_names[pivot_col]

    solution = [0] * num_vars
    for i in range(num_constraints):
        if basis[i] in var_names[:num_vars]:
            idx = var_names.index(basis[i])
            solution[idx] = table[i, -1]

    return table, c, basis, solution

A = [
    [-1, 1],
    [6, 7],
    [2, -3]
]
b = [3, 42, 6]
c = [4, -2]
signs = ['<=', '<=', '<=']
goal = 'min'  # або 'min', залежно від задачі
table, z_row, var_names, basis = build_canonical_table(c, A, b, signs, goal)
print("=== Початкова таблиця ===")
print_table_with_tabulate(table, var_names, basis, z_row)

table, z_row, basis, solution = symbolic_simplex_method_final(table, z_row, var_names, basis, goal)
print("\n=== Результат симплекс-методу ===")
print_table_with_tabulate(Matrix(table), var_names, basis, Matrix(z_row))

sol_str = ", ".join([f"{var_names[i]} = {sp.nsimplify(solution[i], rational=True)}" for i in range(len(var_names))])
F_value = sum(sp.Rational(c[i]) * solution[i] for i in range(len(c)))
F_value = sp.nsimplify(F_value, rational=True)
print(f"Розв’язок: {sol_str}")
print(f"F(x) = {F_value}")

if all(is_integer(val) for val in solution):
    print("\nУсі значення цілі.")
else:
    print("\nПроміжний розв’язок нецілий, починаємо додавати зрізи Гоморі.")

    while True:
        table, z_row, var_names, basis = add_gomory_cut(table, z_row, var_names, basis)

        print("\n=== Симплекс-таблиця після додавання зрізу Гоморі ===")
        print_table_with_tabulate(Matrix(table), var_names, basis, Matrix(z_row))

        table, z_row, basis, solution = dual_simplex(table, z_row, var_names, basis)

        print("\n=== Після розв’язання методом двоїстого симплексу ===")
        print_table_with_tabulate(Matrix(table), var_names, basis, Matrix(z_row))

        if solution is None:
            print("\nНемає допустимого рішення після зрізу Гоморі.")
            break

        if all(is_integer(val) for val in solution):
            print("\nУсі значення цілі.")
            break
        else:
            print("\nПроміжний розв’язок усе ще нецілий, додаємо ще один зріз Гоморі.")

    final_solution = [0] * len(var_names)
    for i in range(len(basis)):
        if var_names.index(basis[i]) < len(var_names):
            final_solution[var_names.index(basis[i])] = table[i, -1]

    sol_str = ", ".join([f"{var_names[i]} = {sp.nsimplify(final_solution[i], rational=True)}" for i in range(len(var_names))])
    F_value = sum(sp.Rational(c[i]) * final_solution[i] for i in range(len(c)))
    F_value = sp.nsimplify(F_value, rational=True)
    print(f"Новий розв’язок: {sol_str}")
    print(f"F(x) = {F_value}")


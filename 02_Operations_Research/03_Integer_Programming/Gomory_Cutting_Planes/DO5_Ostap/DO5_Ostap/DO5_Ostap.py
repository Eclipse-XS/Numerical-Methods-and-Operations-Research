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
            art_var = f"w{len([v for v in var_names if v.startswith('w')]) + 1}"
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

def print_table_with_tabulate(table, var_names, basis, z_row):
    headers = ["Базис"] + [v.upper() for v in var_names] + ["B"]
    data = []
    for i in range(table.rows):
        row = [basis[i]] + [str(nsimplify(cell, rational=True)) for cell in table.row(i)]
        data.append(row)
    z = ["Z"] + [str(nsimplify(e, rational=True)) for e in z_row]
    data.append(z)
    print(tabulate(data, headers=headers, tablefmt="grid", stralign="center"))

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
    num_vars = table.shape[1] - 1  # без вільного члена

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
    # Тільки x1, x2
    target_vars = ["x1", "x2"]
    row_candidates = []

    for i, base_var in enumerate(basis):
        if base_var in target_vars:
            frac_part = get_fractional_part(table[i, -1]) if not is_integer(table[i, -1]) else 0
            row_candidates.append((i, frac_part))

    if not row_candidates:
        return None

    # Знаходимо максимальну дробову частину
    max_row = max(row_candidates, key=lambda x: x[1])

    if max_row[1] == 0:
        return None

    return max_row[0]


def add_gomory_cut(table, c, var_names, basis):
    row_index = find_gomory_row_index(table, var_names, basis)
    if row_index is None:
        return table, c, var_names, basis

    # --- ЛИШЕ x1, x2 ---
    plan_vars = []
    for var in ["x1", "x2"]:
        if var in basis:
            plan_vars.append(var)
        else:
            plan_vars.append("0")
    plan_str = "; ".join(plan_vars)
    print(f"\nОпорний план: ({plan_str})")

    # --- Дробова частина вільного члена ---
    frac_value = nsimplify(get_fractional_part(table[row_index, -1]), rational=True)
    print(f"Max fractional part at {basis[row_index]}: β = {frac_value}")

    # --- Побудова рядка для зрізу ---
    row = table[row_index]
    gomory_coeffs = [-(get_fractional_part(row[i])) for i in range(len(row) - 1)]
    gomory_rhs = -get_fractional_part(row[-1])

    print("\nРахуємо нові коефіцієнти α:")
    alpha_count = 1  # Нумерація α₁, α₂, α₃...
    for idx, val in enumerate(row[:-1]):
        coeff_frac = nsimplify(get_fractional_part(val), rational=True)
        if coeff_frac == 0:
            continue  # Пропускаємо, якщо дробова частина 0
        floor_val = sp.floor(val)
        expr = f"α{alpha_count} = {{ {nsimplify(val, rational=True)} }} = {nsimplify(val, rational=True)} - ⌊{nsimplify(val, rational=True)}⌋ = {nsimplify(val, rational=True)} - ({floor_val}) = {coeff_frac}"
        print(expr)
        alpha_count += 1

    # --- Виведення коефіцієнтів Cut ---
    cut_expr = ""
    for idx, (coeff, var) in enumerate(zip(gomory_coeffs, var_names)):
        if not np.isclose(coeff, 0):
            coeff_frac = nsimplify(coeff, rational=True)
            if cut_expr == "":
                cut_expr += f"{coeff_frac}*{var}"
            else:
                sign = "+" if coeff_frac > 0 else "-"
                cut_expr += f" {sign} {abs(coeff_frac)}*{var}"

    new_var = f"cut{len(var_names)}"
    cut_expr += f" + 1*{new_var}"

    rhs_frac = nsimplify(gomory_rhs, rational=True)

    print("\nCut coefficients:")
    print(cut_expr.strip(), "=", rhs_frac)
    print(f"\nДодано обмеження: {cut_expr.strip()} = {rhs_frac}")

    # --- Формування нової таблиці ---
    gomory_row = [sp.nsimplify(-(get_fractional_part(row[i])), rational=True) for i in range(len(row) - 1)]
    gomory_row += [1, sp.nsimplify(-get_fractional_part(row[-1]), rational=True)]

    new_table = np.hstack((table[:, :-1], np.zeros((table.shape[0], 1)), table[:, -1:]))
    new_table = np.vstack((new_table, np.array(gomory_row, dtype=object)))
    new_c = np.append(c[:-1], [0, c[-1]])
    var_names.append(new_var)
    basis.append(new_var)

    return new_table, new_c, var_names, basis



def print_highlighted_table(table, c, var_names, basis, pivot_row, pivot_col):
    headers = ["Базис"] + [v.upper() for v in var_names] + ["B"]
    data = []

    for i, row in enumerate(table):
        row_data = []
        for j, val in enumerate(row):
            val_str = str(nsimplify(val, rational=True))
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
        val_str = str(nsimplify(val, rational=True))
        if j == pivot_col:
            val_str = colored(val_str, "yellow")
        z_row.append(val_str)
    data.append(["Z"] + z_row)

    print(tabulate(data, headers=headers, tablefmt="grid", stralign="center"))

def dual_simplex(table, c, var_names, basis):
    num_constraints, num_vars = table.shape[0], table.shape[1] - 1
    iter_count = 1

    while True:
        # Якщо всі вільні члени невід'ємні — виходимо
        if all(table[i, -1] >= 0 for i in range(num_constraints)):
            break

        # Опорний рядок — той, де найменший вільний член (b < 0)
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

        # Обираємо стовпець, де мінімальне |c_j/a_j|
        _, pivot_col = min(candidates, key=lambda x: float(x[0].subs(M, 1000) if isinstance(x[0], sp.Basic) else x[0]))

        pivot_element = table[pivot_row, pivot_col]

        print(f"\n=== ІТЕРАЦІЯ {iter_count} ГОМОРІ ===")
        print(f"🔹 Опорний стовпець: {pivot_col + 1}, 🔹 Опорний рядок: {pivot_row + 1}, значення: {nsimplify(pivot_element, rational=True)}")
        print_highlighted_table(table, c, var_names, basis, pivot_row, pivot_col)
        iter_count += 1

        # Нормалізація опорного рядка
        table[pivot_row] = [simplify(x / pivot_element) for x in table[pivot_row]]

        # Обнулення інших рядків
        for i in range(num_constraints):
            if i != pivot_row:
                factor = table[i, pivot_col]
                table[i] = [simplify(table[i][j] - factor * table[pivot_row][j]) for j in range(table.shape[1])]

        # Оновлення оцінок c
        factor = c[pivot_col]
        c = [simplify(c[j] - factor * table[pivot_row][j]) for j in range(len(c))]

        # Оновлення базису
        basis[pivot_row] = var_names[pivot_col]

    # Побудова рішення
    solution = [0] * num_vars
    for i in range(num_constraints):
        if basis[i] in var_names[:num_vars]:
            idx = var_names.index(basis[i])
            solution[idx] = table[i, -1]

    return table, c, basis, solution


# --- ПРИКЛАД ВИКЛИКУ --- #
A = [[1, 2], [2, 4], [6, 1], [1, 3]]
b = [5, 16, 6, 9]
c = [3, 1]
signs = ['>=', '<=', '>=', '<=']
goal = 'min'

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

    # Виводимо фінальний розв'язок
    final_solution = [0] * len(var_names)
    for i in range(len(basis)):
        if var_names.index(basis[i]) < len(var_names):
            final_solution[var_names.index(basis[i])] = table[i, -1]

    sol_str = ", ".join([f"{var_names[i]} = {sp.nsimplify(final_solution[i], rational=True)}" for i in range(len(var_names))])
    F_value = sum(sp.Rational(c[i]) * final_solution[i] for i in range(len(c)))
    F_value = sp.nsimplify(F_value, rational=True)
    print(f"Новий розв’язок: {sol_str}")
    print(f"F(x) = {F_value}")
import numpy as np
import sympy as sp
from sympy import Matrix, simplify
from tabulate import tabulate

M = sp.Symbol("M")
# === Кольори терміналу ===
YELLOW = "\033[43m"
GREEN = "\033[42m"
BLUE = "\033[44m"
RESET = "\033[0m"

# === АВТОМАТИЧНА НОРМАЛІЗАЦІЯ ===
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

# === ФІКС ВІД’ЄМНИХ B В ДВОЇСТІЙ ===
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

# === ДВОЇСТА ЗАДАЧА ===
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

# === ВИВІД ===
def format_expression(coeffs, vars_prefix='x'):
    return " ".join([
        f"{'' if coef < 0 or i == 0 else '+'}{coef}{vars_prefix}{i+1}"
        for i, coef in enumerate(coeffs)
    ])

def format_lp_problem(A, b, c, signs, sense):
    sense_title = "МІНІМІЗАЦІЯ" if sense == 'min' else "МАКСИМІЗАЦІЯ"
    sense_symbol = "→ min" if sense == 'min' else "→ max"
    output = [f"Пряма задача: {sense_title}", f"F = {format_expression(c)} {sense_symbol}"]
    output.append("За умов:")
    for i, row in enumerate(A):
        output.append(f"  {format_expression(row)} {signs[i]} {b[i]}")
    return "\n".join(output)

def format_dual_output(dual):
    sense_symbol = "→ min" if dual['dual_sense'] == 'min' else "→ max"
    output = [f"Двоїста задача:", f"Z = {format_expression(dual['c_dual'], vars_prefix='y')} {sense_symbol}"]
    output.append("За умов:")
    for i in range(len(dual['A_dual'])):
        row_expr = format_expression(dual['A_dual'][i], vars_prefix='y')
        output.append(f"  {row_expr} {dual['constraint_signs'][i]} {dual['b_dual'][i]}")
    return "\n".join(output)

# === SIMPLEX ===
def build_canonical_table(c, A, b, signs, goal='min'):
    A = Matrix(A)
    b = Matrix(b).reshape(len(b), 1)
    var_names = [f"y{i+1}" for i in range(len(c))]
    table = A
    basis = []
    artificial_vars = []

    for i, sign in enumerate(signs):
        if sign == "<=":
            slack = Matrix.zeros(len(signs), 1)
            slack[i] = 1
            table = table.row_join(slack)
            var_names.append(f"y{len(var_names)+1}")
            basis.append(var_names[-1])
        elif sign == ">=":
            surplus = Matrix.zeros(len(signs), 1)
            surplus[i] = -1
            table = table.row_join(surplus)
            var_names.append(f"y{len(var_names)+1}")
            artificial = Matrix.zeros(len(signs), 1)
            artificial[i] = 1
            table = table.row_join(artificial)
            art_name = f"w{len([v for v in var_names if v.startswith('w')]) + 1}"
            var_names.append(art_name)
            basis.append(art_name)
            artificial_vars.append((i, table.cols - 1))
        elif sign == "=":
            artificial = Matrix.zeros(len(signs), 1)
            artificial[i] = 1
            table = table.row_join(artificial)
            art_name = f"w{len([v for v in var_names if v.startswith('w')]) + 1}"
            var_names.append(art_name)
            basis.append(art_name)
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

    data.append(["Z"] + z_data)

    if title:
        print(f"\n📌 {title}")
    print(tabulate(data, headers=headers, tablefmt="grid", stralign="center"))


def symbolic_simplex_method_final(table, c, var_names, basis, goal='min'):
    num_vars = c.shape[0]
    num_constraints = table.shape[0]
    iteration = 1

    print_colored_table(table, var_names, basis, c, title="ПОЧАТКОВА ТАБЛИЦЯ")

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
            print("❌ Задача не має обмеженого розв'язку")
            return None, None, None

        pivot_row = ratios.index(min(ratios))
        pivot_element = table[pivot_row, pivot_col]

        print_colored_table(
            table, var_names, basis, c,
            pivot_row=pivot_row, pivot_col=pivot_col,
            title=f"ІТЕРАЦІЯ {iteration}"
        )

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
        if basis[i] in var_names:
            solution[var_names.index(basis[i])] = table[i, -1]

    return solution, c[-1], basis


def print_table_with_tabulate(table, var_names, basis, z_row):
    headers = ["Базис"] + [v.upper() for v in var_names] + ["B"]
    data = []
    for i in range(table.rows):
        row = [basis[i]] + [str(simplify(cell)) for cell in table.row(i)]
        data.append(row)
    z = ["Z"] + [str(simplify(e)) for e in z_row]
    data.append(z)
    print(tabulate(data, headers=headers, tablefmt="grid", stralign="center"))

# === ВХІДНІ ДАНІ ===
A = [[-1, 1], [6, 7], [2, -3]]
b = [3, 42, 6]
c = [4, -2]
signs = ['<=', '<=', '<=']
sense = 'min'

# === ВИКОНАННЯ ===
print("========== ПРЯМА ЗАДАЧА ==========")
A_norm, b_norm, signs_norm = normalize_primal(A, b, signs, sense)
print(format_lp_problem(A_norm, b_norm, c, signs_norm, sense))

dual = build_dual(A_norm, b_norm, c, sense)

# 🔧 ФІКСИМО ВІД’ЄМНІ B В ДВОЇСТІЙ
dual['A_dual'], dual['b_dual'], dual['constraint_signs'] = fix_negative_rhs_for_dual(
    dual['A_dual'], dual['b_dual'], dual['constraint_signs']
)

print("\n========== ДВОЇСТА ЗАДАЧА ==========")
print(format_dual_output(dual))

print("\n========== РОЗВ'ЯЗАННЯ ДВОЇСТОЇ ==========")
table, z_row, var_names, basis = build_canonical_table(
    dual['c_dual'], dual['A_dual'], dual['b_dual'],
    dual['constraint_signs'], goal=dual['dual_sense']
)
solution, optimal_value, final_basis = symbolic_simplex_method_final(
    table.copy(), z_row.copy(), var_names, basis.copy(), goal=dual['dual_sense']
)

if solution is not None:
    print_table_with_tabulate(table, var_names, final_basis, z_row)
    print("\n========== РЕЗУЛЬТАТ ==========")
    for i, val in enumerate(solution):
        if val != 0:
            print(f"y{i+1} = {simplify(val)}")
    print(f"Оптимальне значення цільової функції: {simplify(optimal_value)}")
else:
    print("🚫 Розв’язок не знайдено або задача не обмежена.")

import numpy as np

M = 99999999

def build_dual(A, b, c, sense='min'):
    A, b, c = np.array(A), np.array(b), np.array(c)
    A_dual, b_dual, c_dual = A.T, c, b
    dual_sense = 'max' if sense == 'min' else 'min'
    constraint_signs = ['<='] * A_dual.shape[0] if sense == 'min' else ['>='] * A_dual.shape[0]
    return {
        'dual_sense': dual_sense,
        'A_dual': A_dual.tolist(),
        'b_dual': b_dual.tolist(),
        'c_dual': c_dual.tolist(),
        'constraint_signs': constraint_signs
    }

def format_expression(coeffs, vars_prefix='x'):
    return " ".join([
        f"{'+' if coef >= 0 and i > 0 else ''}{coef}{vars_prefix}{i+1}"
        for i, coef in enumerate(coeffs)
    ])

def format_lp_problem(A, b, c, sense):
    sense_title = "МІНІМІЗАЦІЯ" if sense == 'min' else "МАКСИМІЗАЦІЯ"
    sense_symbol = "→ min" if sense == 'min' else "→ max"
    output = [f"Пряма задача: {sense_title}", f"F = {format_expression(c)} {sense_symbol}"]
    output.append("За умов:")
    for i, row in enumerate(A):
        output.append(f"  {format_expression(row)} {b[i]}")
    return "\n".join(output)

def format_dual_output(dual):
    sense_symbol = "→ min" if dual['dual_sense'] == 'min' else "→ max"
    output = [f"Двоїста задача:", f"Z = {format_expression(dual['c_dual'], vars_prefix='y')} {sense_symbol}"]
    output.append("За умов:")
    for i in range(len(dual['A_dual'])):
        row_expr = format_expression(dual['A_dual'][i], vars_prefix='y')
        output.append(f"  {row_expr} {dual['constraint_signs'][i]} {dual['b_dual'][i]}")
    return "\n".join(output)

def initialize_simplex_table(goal):
    tables = {
        'max': [
            [2, 1, 3, -1, 0, 1, 4],
            [1, -1, 2, 0, 1, 0, 2]
        ],
        'min': [
            [-2, -1, -3, 1, 0, 0, 4],
            [-1,  1, -2, 0, -1, 1, 2]
        ]
    }
    return np.array(tables[goal], dtype=float)

def build_objective_function(goal):
    return {
        'max': [-2*M + 8, -M + 5, -3*M - 3, M, 0, 0, -4*M],
        'min': [-M - 8, M - 5, -2*M + 3, 0, -M, 0, 2*M]
    }[goal]

def choose_pivot_column(c, goal):
    if goal == 'max':
        return np.argmin(c[:-1]) if np.min(c[:-1]) < 0 else None
    elif goal == 'min':
        positive_indices = [j for j in range(len(c) - 1) if c[j] > 0]
        return max(positive_indices, key=lambda j: c[j]) if positive_indices else None

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

def print_solution(solution, optimal_value, goal, basis):
    if solution is not None:
        prefix = 'y'
        nonzero_vars = [(f"{prefix}{i+1}", solution[i]) for i in range(len(solution)) if solution[i] != 0]
        solution_str = ", ".join(f"{var} = {val:.2f}" for var, val in nonzero_vars) if nonzero_vars else "Всі змінні = 0"
        print("\nОптимальне рішення:")
        print(f"({solution_str})")
        print(f"Оптимальне значення ЦФ: {(optimal_value):.2f}")
    else:
        print("Задача не має розв'язку")

examples = [
    {
        'name': 'Приклад 1',
        'A': [[-2, 1], [-1, -1], [-3, 2]],
        'b': [-8, -5, 3],
        'c': [-4, 2],
        'sense': 'min'
    },
    {
        'name': 'Приклад 2',
        'A': [[2, -1], [1, 1], [3, -2]],
        'b': [8, 5, -3],
        'c': [-4, 2],
        'sense': 'max'
    }
]

for ex in examples:
    print(format_lp_problem(ex['A'], ex['b'], ex['c'], ex['sense']))
    print()
    dual = build_dual(ex['A'], ex['b'], ex['c'], sense=ex['sense'])
    print(format_dual_output(dual))
    goal = dual['dual_sense']
    table = initialize_simplex_table(goal)
    obj_func = build_objective_function(goal)
    solution, optimal_value, basis = simplex_method(table.copy(), obj_func, goal)
    print_solution(solution, optimal_value, goal, basis)
    print("\n" + "="*50 + "\n")

print("Група: ОІ-23")
print("Виконав: Гончаров Вадим Олегович")
print("Номер ЛР: №3")
print("Варіант: №75")
print("=" * 40)
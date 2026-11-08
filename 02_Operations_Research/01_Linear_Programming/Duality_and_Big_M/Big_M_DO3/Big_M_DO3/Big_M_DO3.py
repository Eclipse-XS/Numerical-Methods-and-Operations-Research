import numpy as np

M = 1000

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
    basis = [i + num_vars - num_constraints for i in range(num_constraints)]

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
            return None, None

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

    return solution, c[-1]

def print_solution(solution, optimal_value, goal):
    if solution is not None:
        nonzero_vars = [(f"x{i+1}", solution[i]) for i in range(len(solution)) if solution[i] != 0]
        solution_str = ", ".join(f"{var} = {val:.2f}" for var, val in nonzero_vars) if nonzero_vars else "Всі змінні = 0"
        print("\nОптимальне рішення:")
        print(f"({solution_str})")
        print(f"Оптимальне значення ЦФ: {optimal_value:.2f}")
    else:
        print("Задача не має розв'язку")

for goal in ['max', 'min']:
    table = initialize_simplex_table(goal)
    obj_func = build_objective_function(goal)
    solution, optimal_value = simplex_method(table.copy(), obj_func, goal)
    print(f"\n=== {'Максимізація' if goal == 'max' else 'Мінімізація'} ===")
    print_solution(solution, optimal_value, goal)
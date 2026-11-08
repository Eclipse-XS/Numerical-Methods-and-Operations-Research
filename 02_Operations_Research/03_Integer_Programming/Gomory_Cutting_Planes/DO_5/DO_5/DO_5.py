import numpy as np
from math import fabs

M = 999999999999

def initialize_table():
    return np.array([
        [2, -1, 1, 0, 0, 0, 8],
        [1,  1, 0, 1, 0, 0, 5],
        [-3, 2, 0, 0, -1, 1, 3]
    ], dtype=float)

def choose_pivot_column(c, goal):
    if goal == 'max':
        return np.argmin(c[:-1]) if np.min(c[:-1]) < 0 else None
    return max((j for j in range(len(c) - 1) if c[j] > 0), default=None, key=lambda j: c[j])

def simplex_method(table, c, goal):
    num_vars, num_constraints = len(c) - 1, len(table)
    basis = [col for col in range(num_vars) if list(table[:, col]).count(1) == 1 and list(table[:, col]).count(0) == num_constraints - 1][:num_constraints]

    while True:
        pivot_col = choose_pivot_column(c, goal)
        if pivot_col is None: break
        ratios = [table[i, -1] / table[i, pivot_col] if table[i, pivot_col] > 0 else np.inf for i in range(num_constraints)]
        if all(r == np.inf for r in ratios): return None, None, None, table, c
        pivot_row = np.argmin(ratios)
        basis[pivot_row] = pivot_col
        table[pivot_row] /= table[pivot_row, pivot_col]
        for i in range(num_constraints):
            if i != pivot_row:
                table[i] -= table[i, pivot_col] * table[pivot_row]
        c -= c[pivot_col] * table[pivot_row]

    solution = np.zeros(num_vars)
    for i in range(num_constraints):
        if basis[i] < num_vars:
            solution[basis[i]] = table[i, -1]
    return solution, c[-1], basis, table, c

def is_integer(x): return np.isclose(x, np.round(x))
def get_fractional_part(x): return x - np.floor(x)

def find_gomory_row_index(basis, solution):
    frac_parts = [get_fractional_part(solution[i]) for i in range(2)]
    var = np.argmax(frac_parts)
    return next(((i, frac_parts[var]) for i, b in enumerate(basis) if b == var), (None, None))

def add_gomory_cut(table, c, solution, basis):
    row_index, _ = find_gomory_row_index(basis, solution)
    if row_index is None: return table, c
    row = table[row_index]
    gomory_row = -(row[:-1] - np.floor(row[:-1]))
    gomory_rhs = -(row[-1] - np.floor(row[-1]))
    new_row = np.append(gomory_row, [1, gomory_rhs])
    new_table = np.vstack((np.hstack((table[:, :-1], np.zeros((table.shape[0], 1)), table[:, -1:])), new_row.reshape(1, -1)))
    new_c = np.append(c[:-1], [0, c[-1]])
    return new_table, new_c

def dual_simplex(table, c, basis):
    num_constraints, num_vars = table.shape[0], table.shape[1] - 1
    while True:
        if np.all(table[:, -1] >= 0): break
        pivot_row = np.argmin(table[:, -1])
        row = table[pivot_row, :-1]
        candidates = [(abs(c[j] / row[j]), j) for j in range(num_vars) if row[j] < 0]
        if not candidates: return None, None, None, None
        _, pivot_col = min(candidates)
        table[pivot_row] /= table[pivot_row, pivot_col]
        for i in range(num_constraints):
            if i != pivot_row:
                table[i] -= table[i, pivot_col] * table[pivot_row]
        c -= c[pivot_col] * table[pivot_row]
        basis[pivot_row] = pivot_col
    solution = np.zeros(num_vars)
    for i in range(num_constraints):
        if basis[i] < num_vars:
            solution[basis[i]] = table[i, -1]
    return table, c, basis, solution

def print_solution(solution, optimal_value, title, is_min=False):
    print(f"\n{title}")
    if solution is not None:
        for i in range(2):
            print(f"x{i+1} = {solution[i]:>6.2f}")
        val = -optimal_value if is_min else optimal_value
        print(f"F(x) = {fabs(val):>6.2f}")
    else:
        print("Задача не має розв'язку.")

def print_table(table):
    print("\nНова симплекс-таблиця з обмеженням Гоморі:")
    for row in table:
        print(" | ".join(f"{val:8.2f}" if abs(val) >= 0.005 else "   0.00" for val in row))

print("=" * 40)
print("Група: ОІ-23")
print("Виконав: Гончаров Вадим Олегович")
print("Номер ЛР: №5")
print("Варіант: №75")
print("=" * 40)

table = initialize_table()
c = np.array([-3*M + 4, 2*M - 2, 0, 0, -M, 0, 3*M], dtype=float)

print("\n[ Симплекс-метод ]")
solution, optimal_value, basis, table, c = simplex_method(table.copy(), c.copy(), goal="min")

if solution is None:
    print("Симплекс-метод: задача не має допустимого розв’язку.")
else:
    print_solution(solution, optimal_value, "Опорний план", is_min=True)

    if all(is_integer(solution[i]) for i in range(2)):
        print("Усі значення цілі → оптимальне рішення знайдено.")
    else:
        row_idx, frac_val = find_gomory_row_index(basis, solution)
        print(f"\nМаксимальна дробова частина: {frac_val:.2f}\n→ Додаємо обмеження Гоморі")
        table, c = add_gomory_cut(table, c, solution, basis)
        basis.append(table.shape[1] - 2)
        print_table(table)

        print("\n[ Двоїстий симплекс-метод ]")
        table, c, basis, final_solution = dual_simplex(table, c, basis)
        print_solution(final_solution, -c[-1], "Цілочисловий розв'язок", is_min=True)

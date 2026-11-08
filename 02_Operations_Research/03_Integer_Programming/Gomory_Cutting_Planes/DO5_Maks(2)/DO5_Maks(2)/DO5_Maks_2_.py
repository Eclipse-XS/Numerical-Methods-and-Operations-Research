import numpy as np
from math import floor, fabs
from fractions import Fraction

M = Fraction(999999999999)

def initialize_table():
    return np.array([
        [Fraction(8), Fraction(-5), Fraction(1), 0, 0, 0, 0, 40],
        [Fraction(-6), Fraction(5), 0, 1, 0, 0, 0, 60],
        [Fraction(2), Fraction(1), 0, 0, 1, 0, 0, 14],
        [Fraction(2), Fraction(5), 0, 0, 0, -1, 1, 10]
    ], dtype=object)

def choose_pivot_column(c, goal):
    if goal == 'max':
        return np.argmin([float(ci) for ci in c[:-1]]) if min(c[:-1]) < 0 else None
    return max((j for j in range(len(c) - 1) if c[j] > 0), default=None, key=lambda j: c[j])

def simplex_method(table, c, goal):
    num_vars, num_constraints = len(c) - 1, len(table)
    basis = [col for col in range(num_vars) if list(table[:, col]).count(1) == 1 and list(table[:, col]).count(0) == num_constraints - 1][:num_constraints]

    while True:
        pivot_col = choose_pivot_column(c, goal)
        if pivot_col is None: break
        ratios = [table[i, -1] / table[i, pivot_col] if table[i, pivot_col] > 0 else float('inf') for i in range(num_constraints)]
        if all(r == float('inf') for r in ratios): return None, None, None, table, c
        pivot_row = np.argmin(ratios)
        basis[pivot_row] = pivot_col
        table[pivot_row] = [x / table[pivot_row, pivot_col] for x in table[pivot_row]]
        for i in range(num_constraints):
            if i != pivot_row:
                table[i] = [table[i][j] - table[i][pivot_col] * table[pivot_row][j] for j in range(len(c))]
        c = [c[j] - c[pivot_col] * table[pivot_row][j] for j in range(len(c))]

    solution = [Fraction(0)] * num_vars
    for i in range(num_constraints):
        if basis[i] < num_vars:
            solution[basis[i]] = table[i, -1]
    return solution, c[-1], basis, table, c

def is_integer(x): return x.denominator == 1
def get_fractional_part(x): return x - floor(x)

def find_gomory_row_index(basis, solution):
    frac_parts = [get_fractional_part(solution[i]) for i in range(2)]
    var = np.argmax(frac_parts)
    return next(((i, frac_parts[var]) for i, b in enumerate(basis) if b == var), (None, None))

def add_gomory_cut(table, c, solution, basis):
    row_index, _ = find_gomory_row_index(basis, solution)
    if row_index is None: return table, c
    row = table[row_index]
    gomory_row = [-(x - floor(x)) for x in row[:-1]]
    gomory_rhs = -(row[-1] - floor(row[-1]))
    new_row = gomory_row + [1, gomory_rhs]
    new_table = np.hstack((table[:, :-1], np.zeros((table.shape[0], 1), dtype=object), table[:, -1:]))
    new_table = np.vstack((new_table, np.array(new_row, dtype=object)))
    new_c = c[:-1] + [0, c[-1]]
    return new_table, new_c

def dual_simplex(table, c, basis):
    num_constraints, num_vars = table.shape[0], table.shape[1] - 1
    while True:
        if all(row[-1] >= 0 for row in table): break
        pivot_row = np.argmin([row[-1] for row in table])
        row = table[pivot_row, :-1]
        candidates = [(abs(c[j] / row[j]), j) for j in range(num_vars) if row[j] < 0]
        if not candidates: return None, None, None, None
        _, pivot_col = min(candidates)
        table[pivot_row] = [x / row[pivot_col] for x in table[pivot_row]]
        for i in range(num_constraints):
            if i != pivot_row:
                table[i] = [table[i][j] - table[i][pivot_col] * table[pivot_row][j] for j in range(len(c))]
        c = [c[j] - c[pivot_col] * table[pivot_row][j] for j in range(len(c))]
        basis[pivot_row] = pivot_col
    solution = [Fraction(0)] * num_vars
    for i in range(num_constraints):
        if basis[i] < num_vars:
            solution[basis[i]] = table[i, -1]
    return table, c, basis, solution

def print_solution(solution, optimal_value, title, is_min=False):
    print(f"\n{title}")
    if solution is not None:
        for i in range(2):
            print(f"x{i+1} = {solution[i]}")
        val = -optimal_value if is_min else optimal_value
        print(f"F(x) = {fabs(float(val))}")
    else:
        print("Задача не має розв'язку.")

def print_table(table):
    print("\nНова симплекс-таблиця з обмеженням Гоморі (раціональні дроби):")
    for row in table:
        print(" | ".join(f"{val}" for val in row))

table = initialize_table()
c = [-(2*M + 1), -(5*M + 3), 0, 0, 0, M, 0, -10*M]

print("\n[ Симплекс-метод ]")
solution, optimal_value, basis, table, c = simplex_method(table.copy(), c.copy(), goal="max")

if solution is None:
    print("Симплекс-метод: задача не має допустимого розв’язку.")
else:
    print_solution(solution, optimal_value, "Опорний план", is_min=True)

    if all(is_integer(solution[i]) for i in range(2)):
        print("Усі значення цілі → оптимальне рішення знайдено.")
    else:
        row_idx, frac_val = find_gomory_row_index(basis, solution)
        print(f"\nМаксимальна дробова частина: {frac_val}\n→ Додаємо обмеження Гоморі")
        table, c = add_gomory_cut(table, c, solution, basis)
        basis.append(table.shape[1] - 2)
        print_table(table)

        print("\n[ Двоїстий симплекс-метод ]")
        table, c, basis, final_solution = dual_simplex(table, c, basis)
        print_solution(final_solution, -c[-1], "Цілочисловий розв'язок", is_min=True)

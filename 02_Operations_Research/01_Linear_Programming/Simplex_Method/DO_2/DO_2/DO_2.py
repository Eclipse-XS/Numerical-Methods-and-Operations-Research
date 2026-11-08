import numpy as np

def initialize_simplex_table():
    table = np.array([
        [2, -1, 1, 0, 0, 0, 8],  
        [1, 1, 0, 1, 0, 0, 5],    
        [-3, 2, 0, 0, -1, 1, 3],  
    ], dtype=float)
    
    return table

def simplex_method(table, c):
    num_vars = len(c) - 1  # К-сть змінних
    num_constraints = len(table)  # К-сть обмежень
    basis = [i + num_vars - num_constraints for i in range(num_constraints)]  # Індекси базисних змінних
    
    while True:
        pivot_col = np.argmin(c[:-1])
        if c[pivot_col] >= 0:
            break  # Оптимальне рішення знайдено
        
        ratios = []
        for i in range(num_constraints):
            if table[i, pivot_col] > 0:
                ratios.append(table[i, -1] / table[i, pivot_col])
            else:
                ratios.append(np.inf)
        
        if all(r == np.inf for r in ratios):
            print("Задача не має обмеженого розв’язку.")
            return None, None
        
        pivot_row = np.argmin(ratios)
        basis[pivot_row] = pivot_col  # Оновлення базисної змінної
        
        table[pivot_row, :] /= table[pivot_row, pivot_col]
        
        for i in range(num_constraints):
            if i != pivot_row:
                table[i, :] -= table[i, pivot_col] * table[pivot_row, :]
        
        c -= c[pivot_col] * table[pivot_row, :]

    solution = np.zeros(num_vars)
    for i in range(num_constraints):
        solution[basis[i]] = table[i, -1]  # Встановлення значень базисних змінних
    
    return solution, c[-1]

M = 1000

objective_function_max = [3*M + 4, -2*M - 2, 0, 0, M, 0, -3*M]
objective_function_min = [3*M - 4, -2*M + 2, 0, 0, M, 0, -3*M]

simplex_table = initialize_simplex_table()

solution_max, optimal_value_max = simplex_method(simplex_table.copy(), objective_function_max)
solution_min, optimal_value_min = simplex_method(simplex_table.copy(), objective_function_min)


print("=" * 40)
print("Група: ОІ-23")
print("Виконав: Гончаров Вадим Олегович")
print("Номер ЛР: №2")
print("Варіант: №75")
print("=" * 40)

def print_solution(solution, optimal_value, is_min=False):
    if solution is not None:
        nonzero_vars = [(f"x{i+1}", solution[i]) for i in range(len(solution)) if solution[i] != 0]
        if nonzero_vars:
            solution_str = ", ".join(f"{var} = {val:.2f}" for var, val in nonzero_vars)
        else:
            solution_str = "Всі змінні = 0"

        print("\nОптимальне рішення:")
        print(f"({solution_str})")
        print(f"Оптимальне значення ЦФ: {(-optimal_value if is_min else optimal_value):.2f}")
    else:
        print("Задача не має розв'язку")

print("\n=== Максимізація ===")
print_solution(solution_max, optimal_value_max)

print("\n=== Мінімізація ===")
print_solution(solution_min, optimal_value_min, is_min=True)
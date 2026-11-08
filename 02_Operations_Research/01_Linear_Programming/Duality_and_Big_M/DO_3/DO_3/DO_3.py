import numpy as np

def initialize_table():
    return np.array([
        [2, -1, 1, 0, 0, 8],
        [1,  1, 0, 1, 0, 5],
        [3, -2, 0, 0, 1, -3]
    ], dtype=float)

def dual_simplex(table, c):
    num_constraints, num_columns = table.shape
    num_vars = num_columns - 1
    basis = [2, 3, 4]

    while True:
        rhs = table[:, -1]
        if np.all(rhs >= 0):
            break

        pivot_row = np.argmin(rhs)
        row = table[pivot_row, :-1]

        ratios = []
        for j in range(num_vars):
            if row[j] < 0:
                ratio = abs(c[j] / row[j])
                ratios.append((ratio, j))
        if not ratios:
            return None, None, None, None

        _, pivot_col = min(ratios)

        pivot_element = table[pivot_row, pivot_col]
        table[pivot_row, :] /= pivot_element
        for i in range(num_constraints):
            if i != pivot_row:
                table[i, :] -= table[i, pivot_col] * table[pivot_row, :]

        c -= c[pivot_col] * table[pivot_row, :]
        basis[pivot_row] = pivot_col

    solution = np.zeros(num_vars)
    for i in range(num_constraints):
        solution[basis[i]] = table[i, -1]

    return table, c, basis, solution

def simplex_method(table, c, basis):
    num_constraints, num_columns = table.shape
    num_vars = num_columns - 1

    while True:
        cj = c[:-1]
        if np.all(cj <= 0):
            break

        pivot_col = np.argmax(cj)
        if c[pivot_col] <= 0:
            break

        ratios = []
        for i in range(num_constraints):
            if table[i, pivot_col] > 0:
                ratios.append(table[i, -1] / table[i, pivot_col])
            else:
                ratios.append(np.inf)

        if all(r == np.inf for r in ratios):
            return None, None

        pivot_row = np.argmin(ratios)

        pivot_element = table[pivot_row, pivot_col]
        table[pivot_row, :] /= pivot_element
        for i in range(num_constraints):
            if i != pivot_row:
                table[i, :] -= table[i, pivot_col] * table[pivot_row, :]

        c -= c[pivot_col] * table[pivot_row, :]
        basis[pivot_row] = pivot_col

    solution = np.zeros(num_vars)
    for i in range(num_constraints):
        solution[basis[i]] = table[i, -1]

    return solution, -c[-1]

def print_solution(solution, optimal_value, title, is_min=False):
    print(f"\n{title}")
    if solution is not None:
        for i in range(2):
            print(f"x{i+1} = {solution[i]:>6.2f}")
        val = -optimal_value if is_min else optimal_value
        print(f"F(x) = {val:>6.2f}")
    else:
        print("Задача не має розв'язку.")

if __name__ == '__main__':
    print("=" * 40)
    print(" Група: ОІ-23")
    print(" Виконав: Гончаров Вадим Олегович")
    print(" Номер ЛР: №3")
    print(" Варіант: №75")
    print("=" * 40)

    for mode in ["max", "min"]:
        print("\n" + "=" * 15 + f" ЗАДАЧА {mode.upper()} " + "=" * 15)
        table = initialize_table()
        if mode == "min":
            c = np.array([4, -2, 0, 0, 0, 0], dtype=float)
        else:
            c = np.array([-4, 2, 0, 0, 0, 0], dtype=float)

        print("\n[ Двоїстий симплекс-метод ]")
        result = dual_simplex(table.copy(), c.copy())
        if result[0] is None:
            print("Немає допустимого розв’язку.")
        else:
            table, c, basis, feasible_solution = result
            print_solution(feasible_solution, -c[-1], "Допустимий розв'язок", is_min=(mode=="min"))

            print("\n[ Стандартний симплекс-метод ]")
            solution, optimal_value = simplex_method(table, c, basis)
            print_solution(solution, optimal_value, "Оптимальний розв'язок", is_min=(mode=="min"))
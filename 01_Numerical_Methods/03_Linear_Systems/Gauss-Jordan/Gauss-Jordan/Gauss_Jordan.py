import os
import sympy as sp 

GREEN = '\033[92m'
RESET = '\033[0m'

def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')

def backward_substitution(matrix, variables):
    n = len(matrix)
    final_values = {}
    final_values[variables[-1]] = matrix[-1][-1]
    print(f"x{n} = {final_values[variables[-1]]:.5f}")
    print()

    for i in range(n - 2, -1, -1):
        print(f"Підставляємо знайдені значення в {i+1}-е рівняння:")
        current_value = matrix[i][-1]
        for j in range(i + 1, n):
            current_value -= matrix[i][j] * final_values[variables[j]]

        final_values[variables[i]] = current_value
        print(f"x{i + 1} = {current_value:.5f}")
        print()

    print("\nОстаточні значення змінних:")
    for var in variables:
        print(f"{var} = {final_values[var]:.5f}") 

    return final_values

def print_equations(matrix, variables):
    n = len(matrix)
    for i in range(n):
        equation = ""
        for j in range(n):
            coefficient = matrix[i][j]
            variable = variables[j]
            if j > 0 and coefficient >= 0:
                equation += f" + {coefficient}*{variable}"
            else:
                equation += f" {coefficient}*{variable}"
        
        free_term = matrix[i][-1]
        equation += f" = {free_term}"
        print(equation)
    print()

def print_matrix(matrix, highlighted_row=None):
    print("\nПоточний вигляд матриці:")
    for i in range(len(matrix)):
        formatted_row = [f"{elem:.5f}" for elem in matrix[i]]
        if i == highlighted_row:
            print(GREEN + str(formatted_row) + RESET)
        else:
            print(formatted_row)
    print()

def divide_row(matrix, row, divisor):
    return [element / divisor for element in matrix[row]]

def subtract_rows(matrix, target_row, source_row, multiplier):
    return [matrix[target_row][col] - multiplier * matrix[source_row][col] for col in range(len(matrix[target_row]))]

def forward_elimination(matrix):
    n = len(matrix)
    for i in range(n):
        diagonal_element = matrix[i][i]
        print(f"Нормалізуємо рядок {i + 1}, ділимо на діагональний елемент: {diagonal_element:.5f}")
        matrix[i] = divide_row(matrix, i, diagonal_element)
        print_matrix(matrix, highlighted_row=i)
        input()
        for j in range(i + 1, n):
            factor = matrix[j][i]
            print(f"Обнуляємо елемент в рядку {j + 1} з коефіцієнтом {factor:.5f}")
            matrix[j] = subtract_rows(matrix, j, i, factor)
            print_matrix(matrix, highlighted_row=j)
            input()
    return matrix

def backward_elimination(matrix):
    n = len(matrix)
    for i in range(n - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            factor = matrix[j][i]
            print(f"Обнуляємо елемент в рядку {j + 1} з коефіцієнтом {factor:.5f}")
            matrix[j] = subtract_rows(matrix, j, i, factor)
            print_matrix(matrix, highlighted_row=j) 
            input()
    return matrix

def calculate_absolute_errors(matrix, solutions, variables):
    print("\nОцінка абсолютної похибки:")
    for i, row in enumerate(matrix):
        # Calculate the left-hand side using floating-point substitutions
        calculated_value = sum(float(row[j]) * float(solutions[variables[j]]) for j in range(len(variables)))
        free_term = float(row[-1])
        error = abs(calculated_value - free_term)
        print(f"Рівняння {i + 1}: обчислене значення = {calculated_value:.5f}, вільний член = {free_term}, абсолютна похибка = {error:.5f}")
    print()


def gauss_method(matrix):
    print("=== Метод Гауса ===")
    variables = sp.symbols(f"x1 x2 x3 x4 x5") 
    print_equations(matrix, variables)

    print("=== Прямий хід ===")
    forward_elimination(matrix)
    print("\n=== Зворотній хід ===")
    solutions = backward_substitution(matrix, variables)
    
    calculate_absolute_errors(matrix, solutions, variables)

    return solutions

A = [
    [-1, 22, 18, 22, 15, -8],
    [8, -8, -18, -5, 1, 13],
    [-6, -20, 10, 10, -12, -1],
    [-18, -13, 13, 24, -8, 7],
    [-7, 14, 17, -6, -9, 8]
]

while True:
    clear_screen()
    print("=" * 40)
    print("Група: ОІ-23")
    print("Виконав: Гончаров Вадим Олегович")
    print("Номер ЛР: №1")
    print("Варіант: №4")
    print("=" * 40)
    print("Оберіть метод розв'язання:")
    print("1 - Метод Гауса")
    print("2 - Метод Гауса-Жордана")
    choice = input("Ваш вибір (1 або 2): ")

    if choice in ['1', '2']:
        if choice == '1':
            result = gauss_method(A)
        elif choice == '2':
            result = gauss_jordan_method(A)

            print("\nКінцевий вигляд розширеної матриці:")
            print_matrix(result)

            solutions = [row[-1] for row in result]
            for i, sol in enumerate(solutions):
                print(f"x{i + 1} = {sol:.5f}")

            calculate_absolute_errors(A, {f"x{i + 1}": sol for i, sol in enumerate(solutions)}, sp.symbols(f"x1 x2 x3 x4 x5"))

        break
    else:
        print("Будь ласка, введіть 1 або 2.")
        input("Натисніть Enter, щоб спробувати ще раз...")

print("=" * 40)
print("Група: ОІ-23")
print("Виконав: Гончаров Вадим Олегович")
print("Номер ЛР: №1")
print("Варіант: №4")
print("=" * 40)

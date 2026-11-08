import os
import numpy as np
from tabulate import tabulate

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')

def aitken_scheme_table_compact(x_values, y_values, x):
    n = len(x_values)
    headers = ["xi", "yi", "xi - x"] + [f"L_{','.join([f'i-{j}' for j in range(k, 0, -1)])},i(x)" for k in range(1, n)]
    table = []
    
    L = np.zeros((n, n))
    L[:, 0] = y_values

    for i in range(n):
        row = [x_values[i], y_values[i], x_values[i] - x]
        for j in range(1, n):
            if i >= j:
                L[i, j] = ((x - x_values[i - j]) * L[i, j - 1] - (x - x_values[i]) * L[i - 1, j - 1]) / (x_values[i] - x_values[i - j])
                row.append(L[i, j])
            else:
                row.append("")
        table.append(row)
    
    return headers, table

def experiment_dependence(a, b, x, ns):
    true_value = np.cos(x)
    results = []

    for n in ns:
        h = (b - a) / n
        x_values = np.linspace(a, b, n + 1)
        y_values = np.cos(x_values)

        headers, table = aitken_scheme_table_compact(x_values, y_values, x)
        interpolated_value = table[-1][-1] if table[-1][-1] != "" else table[-2][-1]
        error = abs(interpolated_value - true_value)

        results.append([
            n,
            h,
            interpolated_value,
            true_value,
            error
        ])

    headers = ["Розбиття (n)", "Крок (h)", "Наближене значення", "Точне значення", "Абсолютна Похибка"]
    return headers, results

def main():
    print("=" * 40)
    print("Група: ОІ-23")
    print("Виконав: Гончаров Вадим Олегович")
    print("Номер ЛР: №6")
    print("Варіант: №4")
    print("=" * 40)

    print("Вибір пункту:")
    print("1) Перевірка аналітичного методу")
    print("2) Визначення точності при різній к-сті точок")

    choice = int(input("Обрано: "))
    clear_console()
    if choice == 1:
        x_values = np.arange(0.754, 0.841, 0.0215)
        y_values = np.cos(x_values)
        x = 0.775

        headers, table = aitken_scheme_table_compact(x_values, y_values, x)
        interpolated_value = table[-1][-1] if table[-1][-1] != "" else table[-2][-1]
        true_value = np.cos(x)
        error = abs(interpolated_value - true_value)

        print(f"\nНаближене значення: {interpolated_value:.16f}")
        print(f"Точне значення: {true_value:.16f}")
        print(f"Абсолютна Похибка: {error:.16f}")
        print(tabulate(table, headers=headers, tablefmt="fancy_grid", floatfmt=(".4f", ".16f", ".4f", ".16f", ".16f",".16f",".16f")))

    elif choice == 2:
        a = 0.754
        b = 0.841
        x = 0.775
        ns = [1, 2, 3, 5, 10, 15, 20]

        headers, results = experiment_dependence(a, b, x, ns)
        print(tabulate(results, headers=headers, tablefmt="fancy_grid", floatfmt=(".0f", ".10f", ".16f", ".16f", ".6e")))

if __name__ == "__main__":
    main()

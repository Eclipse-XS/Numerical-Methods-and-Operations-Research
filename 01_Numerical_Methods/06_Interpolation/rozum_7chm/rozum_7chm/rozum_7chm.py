import os
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, simplify, expand, integrate, init_printing
from tabulate import tabulate

init_printing()

x_values = [1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3]
y_values = [1.83, 1.91, 1.99, 2.06, 2.13, 2.2, 2.26, 2.33, 2.39, 2.45, 2.51]

x = symbols('x')

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')

def lagrange_polynomial_with_basis(x_values, y_values):
    n = len(x_values)
    L = 0
    basis_polynomials = []
    for i in range(n):
        l_i = 1
        for j in range(n):
            if i != j:
                l_i *= (x - x_values[j]) / (x_values[i] - x_values[j])
        l_i = expand(simplify(l_i))
        basis_polynomials.append(l_i)
        L += l_i * y_values[i]
    L = expand(simplify(L))
    return L, basis_polynomials

polynomial, basis_polynomials = lagrange_polynomial_with_basis(x_values, y_values)

a, b = 1, 3
integral_high_precision = float(integrate(polynomial, (x, a, b)))

def calculate_errors(exact, approx):
    absolute_error = abs(exact - approx)
    relative_error = (absolute_error / exact) * 100
    return absolute_error, relative_error

def method_of_rectangles(f, a, b, eps, method='left'):
    eps=10**-3
    n = 1
    integral_old = 0
    while True:
        h = (b - a) / n
        if method == 'left':
            x_values = np.linspace(a, b - h, n)
        elif method == 'right':
            x_values = np.linspace(a + h, b, n)
        elif method == 'middle':
            x_values = np.linspace(a + h / 2, b - h / 2, n)
        else:
            raise ValueError("Невідомий тип методу прямокутників")
        integral = h * np.sum(f(x_values))
        if abs(integral - integral_old) < eps:
            return integral, h, n
        integral_old = integral
        n *= 2
        if n > 1e6:
            raise RuntimeError("Метод не збігається")

def trapezoidal_method(f, a, b, eps):
    n = 1
    integral_old = 0
    while True:
        h = (b - a) / n
        x_values = np.linspace(a, b, n + 1)
        y_values = f(x_values)
        integral = h * (0.5 * y_values[0] + np.sum(y_values[1:-1]) + 0.5 * y_values[-1])
        if abs(integral - integral_old) < eps:
            return integral, h, n
        integral_old = integral
        n *= 2
        if n > 1e6:
            raise RuntimeError("Метод не збігається")

def simpson_method(f, a, b, eps):
    n = 2
    integral_old = 0
    while True:
        h = (b - a) / n
        x_values = np.linspace(a, b, n + 1)
        y_values = f(x_values)
        integral = h / 3 * (y_values[0] + 4 * np.sum(y_values[1:-1:2]) + 2 * np.sum(y_values[2:-2:2]) + y_values[-1])
        if abs(integral - integral_old) < eps:
            return integral, h, n
        integral_old = integral
        n *= 2
        if n > 1e6:
            raise RuntimeError("Метод не збігається")

def three_eighths_method(f, a, b, eps):
    n = 3
    integral_old = 0
    while True:
        h = (b - a) / n
        x_values = np.linspace(a, b, n + 1)
        y_values = f(x_values)
        integral = 3 * h / 8 * (y_values[0] + 3 * np.sum(y_values[1:-1:3]) + 3 * np.sum(y_values[2:-1:3]) + 2 * np.sum(y_values[3:-3:3]) + y_values[-1])
        if abs(integral - integral_old) < eps:
            return integral, h, n
        integral_old = integral
        n += 3
        if n > 1e6:
            raise RuntimeError("Метод не збігається")

def polynomial_function(x_values):
    return np.array([float(polynomial.subs(x, val)) for val in x_values])

def display_results_as_table(method, name, method_type=None):
    eps = 10**-6
    results = []
    try:
        if method_type:
            result, h, n = method(polynomial_function, a, b, method=method_type, eps=eps)
        else:  
            result, h, n = method(polynomial_function, a, b, eps=eps)
        abs_error, rel_error = calculate_errors(integral_high_precision, result)
        results.append([eps, h, n, result, abs_error, rel_error])
    except RuntimeError as e:
        results.append([eps, "-", "-", "Не збігся", "-", "-"])
    
    headers = ["Точність", "Крок h", "Розбиття n", "Результат", "Абсолютна похибка", "Відносна похибка (%)"]
    print(f"\nТаблиця результатів для методу: {name}")
    print(tabulate(results, headers=headers, tablefmt="pipe", floatfmt=(".6f", ".8f", ".0f", ".16f", ".10f", ".10f")))

while True:
    clear_console()
    print("--- МЕНЮ ---")
    print("1: Переглянути базисні поліноми та поліном Лагранжа.")
    print("2: Обчислити значення функції для заданого x (з інтервалу [1, 3]).")
    print("3: Побудувати графік полінома Лагранжа.")
    print("4: Обчислити інтеграл (вибір методу).")
    print("5: Вийти.")
    
    choice = input("Введіть номер опції: ")

    if choice == "1":
        clear_console()
        print("\n" + "=" * 70)
        print("                БАЗИСНІ ПОЛІНОМИ l_i(x)")
        print("=" * 70)
        for i, l_i in enumerate(basis_polynomials):
            print(f"\nБазисний поліном l_{i}(x):")
            print("-" * 70)
            print(l_i)
            print("-" * 70)
        print("\n" + "=" * 70)
        print("               ПОЛІНОМ ЛАГРАНЖА L(x)")
        print("=" * 70)
        print(polynomial)
        print("=" * 70)
        input("\nНатисніть Enter, щоб повернутися до меню...")

    elif choice == "2":
        clear_console()
        try:
            user_input = input("Введіть значення x (дійсне число від 1 до 3): ")
            x_value = float(user_input)
            if x_value < 1 or x_value > 3:
                print("Значення x має належати інтервалу [1, 3]. Спробуйте ще раз.")
                continue
            result = polynomial.subs(x, x_value)
            print(f"\nL({x_value}) = {result}")
        except ValueError:
            print("Некоректне значення. Введіть дійсне число.")
        input("\nНатисніть Enter, щоб повернутися до меню...")

    elif choice == "3":
        clear_console()
        x_plot = np.linspace(1, 3, 500)
        y_plot = [polynomial.subs(x, val) for val in x_plot]

        plt.figure(figsize=(10, 6))
        plt.plot(x_plot, y_plot, label="Поліном Лагранжа", color="blue")
        plt.scatter(x_values, y_values, color="red", label="Точки даних")
        plt.title("Інтерполяційний поліном Лагранжа")
        plt.xlabel("x")
        plt.ylabel("L(x)")
        plt.legend()
        plt.grid()
        plt.show()
        input("\nНатисніть Enter, щоб повернутися до меню...")

    elif choice == "4":
        clear_console()
        print("--- ВИБІР МЕТОДУ ОБЧИСЛЕННЯ ІНТЕГРАЛА ---")
        print("1: Точне значення (аналітичне).")
        print("2: Метод прямокутників (лівий, правий, середній).")
        print("3: Метод трапецій.")
        print("4: Метод Сімпсона.")
        print("5: Метод трьох восьмих.")
        method_choice = input("Введіть номер методу: ")

        if method_choice == "1":
            clear_console()
            print(f"Точне значення інтеграла: {integral_high_precision:.16f}")
        elif method_choice == "2":
            clear_console()
            print("--- ВИБІР ПІДМЕТОДУ ПРЯМОКУТНИКІВ ---")
            print("1: Ліві прямокутники.")
            print("2: Праві прямокутники.")
            print("3: Середні прямокутники.")
            sub_method_choice = input("Введіть номер підметоду: ")
            if sub_method_choice == "1":
                display_results_as_table(method_of_rectangles, "Ліві прямокутники", method_type='left')
            elif sub_method_choice == "2":
                display_results_as_table(method_of_rectangles, "Праві прямокутники", method_type='right')
            elif sub_method_choice == "3":
                display_results_as_table(method_of_rectangles, "Середні прямокутники", method_type='middle')
            else:
                print("Некоректний вибір підметоду.")
        elif method_choice == "3":
            clear_console()
            display_results_as_table(trapezoidal_method, "Метод трапецій")
        elif method_choice == "4":
            clear_console()
            display_results_as_table(simpson_method, "Метод Сімпсона")
        elif method_choice == "5":
            clear_console()
            display_results_as_table(three_eighths_method, "Метод трьох восьмих")
        else:
            print("Некоректний вибір методу.")
        input("\nНатисніть Enter, щоб повернутися до меню...")

    elif choice == "5":
        print("Вихід із програми.")
        break

    else:
        print("Неправильний вибір. Спробуйте ще раз.")

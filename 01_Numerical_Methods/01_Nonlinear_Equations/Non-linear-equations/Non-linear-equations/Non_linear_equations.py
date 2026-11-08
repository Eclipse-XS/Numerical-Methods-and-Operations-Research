import os
import sympy as sp

x = sp.Symbol('x')
f = 2 * sp.log(x) - 1 / x
f_der = sp.diff(f, x)
f_der2 = sp.diff(f_der, x)
g = sp.sqrt(sp.exp(1 / x))
g_der = sp.diff(g, x)

def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')

def validate_sign_change(func, a, b):
    if func.subs(x, a) * func.subs(x, b) >= 0:
        raise ValueError("Метод не можна застосувати: f(a) * f(b) ≥ 0.")

def bisection(a, b, acc):
    validate_sign_change(f, a, b)
    iterations = 0
    while (b - a) / 2 > acc:
        mid = (a + b) / 2
        if f.subs(x, mid) == 0:
            return mid, iterations
        elif f.subs(x, a) * f.subs(x, mid) < 0:
            b = mid
        else:
            a = mid
        iterations += 1
    return (a + b) / 2, iterations

def newton(x0, acc):
    if f.subs(x, x0) * f_der2.subs(x, x0) <= 0:
        raise ValueError("Метод Ньютона не можна застосувати: не виконується f(x0)*f''(x0) > 0.")
    validate_sign_change(f, x0 - 1, x0 + 1)
    iterations = 0
    x_prev = x0
    x_new = x0
    while True:
        fx = f.subs(x, x_new)
        fpx = f_der.subs(x, x_new)
        if fpx == 0:
            raise ValueError("Метод Ньютона не можна застосувати: похідна дорівнює 0.")
        x_prev = x_new
        x_new = x_new - fx / fpx
        iterations += 1
        if abs(x_new - x_prev) < acc:
            break
    return x_new, iterations

def chord(x0, x1, acc):
    validate_sign_change(f, x0, x1)
    if f_der2 is None or not f_der2.has(x):
        raise ValueError("Метод хорд не можна застосувати: друга похідна f''(x) не існує.")
    iterations = 0
    f_x0 = f.subs(x, x0)
    while True:
        f_x1 = f.subs(x, x1)
        x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0)
        if abs(x2 - x1) < acc:
            return x2, iterations + 1
        x1 = x2
        iterations += 1


def simple_iteration(g, x0, acc):
    validate_sign_change(f, x0 - 1, x0 + 1)
    if abs(g_der.subs(x, x0)) >= 1:
        raise ValueError("Метод простої ітерації не можна застосувати: |g'(x0)| ≥ 1.")
    iterations = 0
    while True:
        x1 = g.subs(x, x0)
        if abs(x1 - x0) < acc:
            return x1, iterations + 1
        x0 = x1
        iterations += 1

def choose_accuracy():
    print("Оберіть точність:")
    for i in range(1, 11):
        acc = 10 ** -i
        print(f"{i} - {acc:.10f}")
    choice = int(input("Введіть номер точності: "))
    if 1 <= choice <= 10:
        return 10 ** -choice
    else:
        print("Неправильний вибір точності.")
        return choose_accuracy()

def choose_method(tol):
    print("=" * 40)
    print("Група: ОІ-23")
    print("Виконав: Гончаров Вадим Олегович")
    print("Номер ЛР: №1")
    print("Варіант: №4")
    print("=" * 40)

    print(f"Поточна точність: {tol:.10f}")
    print("Оберіть метод для розв'язання рівняння:")
    print("1 - Метод половинного ділення")
    print("2 - Метод Ньютона")
    print("3 - Метод хорд")
    print("4 - Метод простої ітерації")
    print("5 - Обрати точність")
    print("6 - Проаналізувати всі методи на різних рівнях точності")
    
    choice = int(input("Введіть номер методу: "))
    try:
        if choice == 1:
            result, iterations = bisection(1, 2, tol)
            print(f"Розв'язок (половинне ділення): x = {result}, ітерацій: {iterations}")
        elif choice == 2:
            result, iterations = newton(1.0, tol)
            print(f"Розв'язок (метод Ньютона): x = {result}, ітерацій: {iterations}")
        elif choice == 3:
            result, iterations = chord(1.0, 2.0, tol)
            print(f"Розв'язок (метод хорд): x = {result}, ітерацій: {iterations}")
        elif choice == 4:
            result, iterations = simple_iteration(g, 1.5, tol)
            print(f"Розв'язок (проста ітерація): x = {result}, ітерацій: {iterations}")
        elif choice == 5:
            clear_screen()
            new_acc = choose_accuracy()
            clear_screen()
            choose_method(new_acc)
        elif choice == 6:
            analyze_iterations_for_all_methods()
        else:
            print("Неправильний вибір методу.")
            choose_method(tol)
    except ValueError as e:
        print(f"Помилка: {e}")
        choose_method(tol)

def analyze_iterations_for_all_methods():
    accuracies = [10**-i for i in range(1, 11)]
    methods = {
        "Метод половинного ділення": bisection,
        "Метод Ньютона": newton,
        "Метод хорд": chord,
        "Проста ітерація": simple_iteration
    }
    for acc in accuracies:
        print(f"\nПоточна точність: {acc:.10f}")
        for method_name, method in methods.items():
            try:
                if method_name == "Метод половинного ділення":
                    result, iterations = method(1, 2, acc)
                elif method_name == "Метод Ньютона":
                    result, iterations = method(1.0, acc)
                elif method_name == "Метод хорд":
                    result, iterations = method(1.0, 2.0, acc)
                elif method_name == "Проста ітерація":
                    result, iterations = method(g, 1.5, acc)
                print(f"{method_name}: x = {result}, ітерацій: {iterations}")
            except Exception as e:
                print(f"{method_name}: не вдалося виконати через {str(e)}")
    print()

accuracy = choose_accuracy()
clear_screen()
while True:
    choose_method(accuracy)

import sympy as sp
from tabulate import tabulate
import numpy as np

def newton_method_for_two_equations(x_0, y_0, tol, max_iter=100):
    iter_count = 0

    J = sp.Matrix([
        [sp.diff(f, x), sp.diff(f, y)],
        [sp.diff(g, x), sp.diff(g, y)]
    ])

    x_k = x_0
    y_k = y_0

    while iter_count < max_iter:
        f_val = f.subs({x: x_k, y: y_k}).evalf()
        g_val = g.subs({x: x_k, y: y_k}).evalf()

        J_x_y = J.subs({x: x_k, y: y_k}).evalf()
        det_J = J_x_y.det()

        if det_J == 0:
            print("Помилка: визначник матриці Якобі дорівнює нулю. Неможливо продовжити.")
            return None, None, None

        J_inv = 1 / det_J
        x_new = x_k - J_inv * (J_x_y[1, 1] * f_val - J_x_y[0, 1] * g_val)
        y_new = y_k - J_inv * (J_x_y[0, 0] * g_val - J_x_y[1, 0] * f_val)

        if abs(x_new - x_k) < tol and abs(y_new - y_k) < tol:
            print(f"\nЗбіжність досягнута за {iter_count + 1} ітерацій з точністю {tol}")
            return x_new, y_new, iter_count + 1

        x_k, y_k = x_new, y_new
        iter_count += 1

    print("Збіжність не досягнута за максимальну кількість ітерацій.")
    return None, None, None


def simple_iteration_method_for_two_equations(x_0, y_0, tol, max_iter=100):
    phi_1 = (1.6 - sp.sin(y)) / 2
    phi_2 = 1 - sp.cos(x - 1)

    x_k = x_0
    y_k = y_0

    iter_count = 0

    while iter_count < max_iter:
        x_new = phi_1.subs({x: x_k, y: y_k}).evalf()
        y_new = phi_2.subs({x: x_k, y: y_k}).evalf()

        print(f"Ітерація {iter_count + 1}: x = {x_new}, y = {y_new}")

        if abs(x_new - x_k) < tol and abs(y_new - y_k) < tol:
            print(f"\nЗбіжність досягнута за {iter_count + 1} ітерацій з точністю {tol}")
            return x_new, y_new, iter_count + 1

        x_k = x_new
        y_k = y_new
        iter_count += 1

    print("Збіжність не досягнута за максимальну кількість ітерацій.")
    return None, None, None


def research_newton_method_for_two_equations(x_k, y_k, tols):
    data = []

    for tol in tols:
        x_k1, y_k1, iter_count = newton_method_for_two_equations(x_k, y_k, tol)
        data.append([tol, iter_count, x_k1, y_k1])

    print("\nРезультати для методу Ньютона:")
    print(tabulate(data, headers=["Точність", "Ітерації", "x", "y"], tablefmt="grid"))


def research_simple_iteration_method_for_two_equations(x_0, y_0, tols):
    data = []

    for tol in tols:
        x_k, y_k, iter_count = simple_iteration_method_for_two_equations(x_0, y_0, tol)
        data.append([tol, iter_count, x_k, y_k])

    print("\nРезультати для методу простої ітерації:")
    print(tabulate(data, headers=["Точність", "Ітерації", "x", "y"], tablefmt="grid"))


if __name__ == '__main__':
    x, y = sp.symbols('x y')
    delta_x, delta_y = sp.symbols('delta_x delta_y')

    f = sp.cos(x - 1) + y - 1
    g = sp.sin(y) + 2 * x - 1.6

    x_0 = 0.8
    y_0 = 0

    tols = [10**(-x) for x in range(1, 11)]

    print("Оберіть метод для застосування:")
    print("1 - Метод Ньютона")
    print("2 - Метод простої ітерації")
    choice = input("Введіть номер методу, який хочете використати: ").strip()

    if choice == "1":
        research_newton_method_for_two_equations(x_0, y_0, tols)
    elif choice == "2":
        research_simple_iteration_method_for_two_equations(x_0, y_0, tols)
    else:
        print("Неправильний вибір. Будь ласка, перезапустіть програму і оберіть коректний варіант.")

print("=" * 40)
print("Група: ОІ-23")
print("Виконав: Гончаров Вадим Олегович")
print("Номер ЛР: №5")
print("Варіант: №4")
print("=" * 40)

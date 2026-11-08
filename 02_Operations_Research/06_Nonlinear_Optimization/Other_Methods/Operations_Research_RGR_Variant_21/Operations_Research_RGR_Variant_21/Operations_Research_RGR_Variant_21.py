import sympy as sp

x = sp.Symbol('x')
f_expr = 80 * x - 30 * x**2 - 0.25 * x**4
f_lambdified = sp.lambdify(x, f_expr)

x1 = 1.0
dx = 0.5
x2 = x1 + dx
f1 = f_lambdified(x1)
f2 = f_lambdified(x2)

x3 = x1 + 2 * dx if f2 > f1 else x1 - dx
f3 = f_lambdified(x3)

points = [(x1, f1), (x2, f2), (x3, f3)]
points.sort()

epsilon = 0.01
iteration = 1

def show_substitution(val):
    val_str = f"{val:.5f}"
    sub_expr = f"80·{val_str} - 30·({val_str})² - 0.25·({val_str})⁴"
    result = f_lambdified(val)
    return f"f({val_str}) = {sub_expr} = {result:.5f}"

while True:
    (x1, f1), (x2, f2), (x3, f3) = points

    print(f"Ітерація {iteration}:\n")

    print("Обчислення значень функції:")
    print(show_substitution(x1))
    print(show_substitution(x2))
    print(show_substitution(x3))

    f_values = [(x1, f1), (x2, f2), (x3, f3)]
    x_max, f_max = max(f_values, key=lambda t: t[1])
    print(f"\nF_max = f({x_max:.5f}) = {f_max:.5f} (найбільше з трьох значень)\n")

    num = (x2**2 - x3**2)*f1 + (x3**2 - x1**2)*f2 + (x1**2 - x2**2)*f3
    den = (x2 - x3)*f1 + (x3 - x1)*f2 + (x1 - x2)*f3

    print("Обчислення чисельника та знаменника:")
    print(f"Чисельник = ({x2**2:.5f} - {x3**2:.5f})*{f1:.5f} + ({x3**2:.5f} - {x1**2:.5f})*{f2:.5f} + "
          f"({x1**2:.5f} - {x2**2:.5f})*{f3:.5f} = {num:.5f}")
    print(f"Знаменник = ({x2:.5f} - {x3:.5f})*{f1:.5f} + ({x3:.5f} - {x1:.5f})*{f2:.5f} + "
          f"({x1:.5f} - {x2:.5f})*{f3:.5f} = {den:.5f}")

    if den == 0:
        print("Ділення на нуль. Завершення.")
        break

    x_star = 0.5 * num / den
    f_star = f_lambdified(x_star)

    print(f"\nx* = 0.5 * ({num:.5f} / {den:.5f}) = {x_star:.5f}")
    print(show_substitution(x_star))

    print("\nПеревірка збіжності:")
    diff = abs(x_star - x_max)
    print(f"|{x_star:.5f} - {x_max:.5f}| = {diff:.5f} < {epsilon} → {'Так' if diff < epsilon else 'Ні'}")

    if diff < epsilon:
        print("\nКритерій збіжності виконано. Завершення.\n")
        break
    else:
        print("\nКритерій збіжності не виконано. Переходимо до наступної ітерації.\n")

    new_points = points + [(x_star, f_star)]
    new_points.sort()
    best_index = [pt[0] for pt in new_points].index(x_star)

    if best_index == 0:
        points = new_points[:3]
    elif best_index == len(new_points) - 1:
        points = new_points[-3:]
    else:
        points = new_points[best_index - 1: best_index + 2]

    iteration += 1

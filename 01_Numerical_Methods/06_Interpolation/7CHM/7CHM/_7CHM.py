import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from tabulate import tabulate

x = np.array([1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3])
y = np.array([0.19, 0.35, 0.52, 0.7, 0.89, 1.08, 1.28, 1.48, 1.67, 1.87, 2.07])

cs = CubicSpline(x, y, bc_type='natural')
x_dense = np.linspace(min(x), max(x), 500)  
y_dense = cs(x_dense)


a, b = 1, 3
integral_high_precision = cs.integrate(a, b)

print(f"Точне значення інтеграла (з точністю до 10^-16): {integral_high_precision:.16f}")


def calculate_errors(exact, approx):
    absolute_error = abs(exact - approx)
    relative_error = (absolute_error / exact) * 100  
    return absolute_error, relative_error

def method_of_rectangles(func, a, b, method='left', eps=1e-6):
    n = 1
    while True:
        h = (b - a) / n
        if method == 'left':
            x_points = np.linspace(a, b - h, n)
            integral = h * np.sum(func(x_points))
        elif method == 'right':
            x_points = np.linspace(a + h, b, n)
            integral = h * np.sum(func(x_points))
        elif method == 'middle':
            x_points = np.linspace(a + h / 2, b - h / 2, n)
            integral = h * np.sum(func(x_points))
        else:
            raise ValueError("Invalid method! Choose 'left', 'right', or 'middle'.")
        
        n_new = 2 * n
        h_new = (b - a) / n_new
        if method == 'left':
            x_points_new = np.linspace(a, b - h_new, n_new)
            integral_new = h_new * np.sum(func(x_points_new))
        elif method == 'right':
            x_points_new = np.linspace(a + h_new, b, n_new)
            integral_new = h_new * np.sum(func(x_points_new))
        elif method == 'middle':
            x_points_new = np.linspace(a + h_new / 2, b - h_new / 2, n_new)
            integral_new = h_new * np.sum(func(x_points_new))
        
        if abs(integral_new - integral) < eps:
            return integral, h, n
        n = n_new

def display_rectangles_results_as_table():
    print(f"\nТочне значення інтеграла: {integral_high_precision:.16f}")
    a, b = 1, 3
    eps_values = [10**-i for i in range(1, 7)]
    
    results = {
        "Ліві прямокутники": [],
        "Праві прямокутники": [],
        "Середні прямокутники": []
    }
    
    for eps in eps_values:
        left_result, left_h, left_n = method_of_rectangles(cs, a, b, method='left', eps=eps)
        left_abs_error, left_rel_error = calculate_errors(integral_high_precision, left_result)
        results["Ліві прямокутники"].append([
            f"{eps:.6f}",
            f"{left_h:.16f}",
            left_n,
            f"{left_result:.16f}",
            f"{left_abs_error:.16f}",
            f"{left_rel_error:.10f}"
        ])
        
        right_result, right_h, right_n = method_of_rectangles(cs, a, b, method='right', eps=eps)
        right_abs_error, right_rel_error = calculate_errors(integral_high_precision, right_result)
        results["Праві прямокутники"].append([
            f"{eps:.6f}",
            f"{right_h:.16f}",
            right_n,
            f"{right_result:.16f}",
            f"{right_abs_error:.16f}",
            f"{right_rel_error:.10f}"
        ])
        
        middle_result, middle_h, middle_n = method_of_rectangles(cs, a, b, method='middle', eps=eps)
        middle_abs_error, middle_rel_error = calculate_errors(integral_high_precision, middle_result)
        results["Середні прямокутники"].append([
            f"{eps:.6f}",
            f"{middle_h:.16f}",
            middle_n,
            f"{middle_result:.16f}",
            f"{middle_abs_error:.16f}",
            f"{middle_rel_error:.10f}"
        ])
    
    for method_name, method_results in results.items():
        print(f"\nТаблиця результатів для {method_name}:")
        table = tabulate(
            method_results,
            headers=["Точність", "Крок h", "Розбиття n", "Результат", "Абсолютна похибка", "Відносна похибка (%)"],
            tablefmt='fancy_grid',
            floatfmt=(".6f", ".12f", ".2f", ".16f", ".16f",".16f",".16f")
        )
        print(table)

def trapezoidal_method(func, a, b, eps=1e-6):
    n = 1
    while True:
        h = (b - a) / n
        x_points = np.linspace(a, b, n + 1)
        y_points = func(x_points)
        integral = (h / 2) * (y_points[0] + 2 * np.sum(y_points[1:-1]) + y_points[-1])
        
        n_new = 2 * n
        h_new = (b - a) / n_new
        x_points_new = np.linspace(a, b, n_new + 1)
        y_points_new = func(x_points_new)
        integral_new = (h_new / 2) * (y_points_new[0] + 2 * np.sum(y_points_new[1:-1]) + y_points_new[-1])
        
        if abs(integral_new - integral) < eps:
            return integral_new, h, n
        n = n_new

def simpson_method(func, a, b, eps=1e-6):
    n = 2  
    while True:
        h = (b - a) / n
        x_points = np.linspace(a, b, n + 1)
        y_points = func(x_points)
        integral = (h / 3) * (
            y_points[0] +
            4 * np.sum(y_points[1:-1:2]) +
            2 * np.sum(y_points[2:-2:2]) +
            y_points[-1]
        )
        
        n_new = 2 * n
        h_new = (b - a) / n_new
        x_points_new = np.linspace(a, b, n_new + 1)
        y_points_new = func(x_points_new)
        integral_new = (h_new / 3) * (
            y_points_new[0] +
            4 * np.sum(y_points_new[1:-1:2]) +
            2 * np.sum(y_points_new[2:-2:2]) +
            y_points_new[-1]
        )
        
        if abs(integral_new - integral) < eps:
            return integral_new, h, n
        n = n_new

def three_eighths_method(func, a, b, eps=1e-6):
    n = 3  
    while True:
        h = (b - a) / n
        x_points = np.linspace(a, b, n + 1)
        y_points = func(x_points)
        integral = (3 * h / 8) * (
            y_points[0] +
            3 * np.sum(y_points[1:n:3]) +
            3 * np.sum(y_points[2:n:3]) +
            2 * np.sum(y_points[3:n-1:3]) +
            y_points[-1]
        )
        
        n_new = n + 3
        h_new = (b - a) / n_new
        x_points_new = np.linspace(a, b, n_new + 1)
        y_points_new = func(x_points_new)
        integral_new = (3 * h_new / 8) * (
            y_points_new[0] +
            3 * np.sum(y_points_new[1:n_new:3]) +
            3 * np.sum(y_points_new[2:n_new:3]) +
            2 * np.sum(y_points_new[3:n_new-1:3]) +
            y_points_new[-1]
        )
        
        if abs(integral_new - integral) < eps:
            return integral_new, h, n
        n = n_new

def display_results_as_table(method, name, method_type=None):
    print(f"\nТочне значення інтеграла: {integral_high_precision:.16f}")
    eps_values = [10**-i for i in range(1, 7)]

    results = []
    for eps in eps_values:
        if method_type:
            result, h, n = method(cs, a, b, method=method_type, eps=eps)
        else:  
            result, h, n = method(cs, a, b, eps=eps)
        abs_error, rel_error = calculate_errors(integral_high_precision, result)
        results.append([
            eps,
            h,
            n,
            result,
            abs_error,
            rel_error
        ])

    headers = ["Точність", "Крок h", "Розбиття n", "Результат", "Абсолютна похибка", "Відносна похибка (%)"]
    print(tabulate(
        results,
        headers=headers,
        tablefmt="fancy_grid",
        floatfmt=(".6f", ".8f", ".0f", ".16f", ".10f", ".10f")
    ))

print("=" * 40)
print("Група: ОІ-23")
print("Виконав: Гончаров Вадим Олегович")
print("Номер ЛР: №7")
print("Варіант: №4")
print("=" * 40)

print("Оберіть дію:")
print("1) Метод прямокутників")
print("2) Метод трапецій")
print("3) Метод Сімпсона")
print("4) Метод трьох восьмих")
print("5) Вивести графік кубічного сплайна")
print("6) Показати рівняння кубічного сплайна")
choice = int(input("Ваш вибір: "))

if choice == 1:
    print("\nТаблиця результатів для методу прямокутників:")
    display_rectangles_results_as_table()
elif choice == 2:
    print("\nТаблиця результатів для методу трапецій:")
    display_results_as_table(trapezoidal_method, 'Трапеції')
elif choice == 3:
    print("\nТаблиця результатів для методу Сімпсона:")
    display_results_as_table(simpson_method, 'Сімпсон')
elif choice == 4:
    print("\nТаблиця результатів для методу трьох восьмих:")
    display_results_as_table(three_eighths_method, 'Три восьмі')
elif choice == 5:
    plt.figure(figsize=(10, 6))
    plt.plot(x_dense, y_dense, label='Кубічний сплайн', color='black', linewidth=2)  
    plt.scatter(x, y, color='red', label='Задані точки (x, y)', zorder=5)  
    plt.title("Графік кубічного сплайну")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.grid(True)
    plt.show()
elif choice == 6:
    print("\nРівняння кубічного сплайну для кожного інтервалу:")
    for i in range(len(cs.c.T)):
        c3, c2, c1, c0 = cs.c.T[i] 
        print(f"[{x[i]:.1f}, {x[i+1]:.1f}]: S(x) = "
              f"({c3:.6f}) * (x - {x[i]:.1f})^3 + "
              f"({c2:.6f}) * (x - {x[i]:.1f})^2 + "
              f"({c1:.6f}) * (x - {x[i]:.1f}) + "
              f"({c0:.6f})")


else:
    print("Неправильний вибір!")

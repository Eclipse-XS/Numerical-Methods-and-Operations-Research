import numpy as np
import matplotlib.pyplot as plt

# Створення координат для осей X та Y
x = np.linspace(0, 10, 400)
y = np.linspace(0, 10, 400)

# Створення сітки точок
X, Y = np.meshgrid(x, y)

# Перевірка умов для кожної точки
condition1 = X - Y >= -3  # для x - y >= -3
condition2 = 6 * X + 7 * Y <= 42  # для 6x + 7y <= 42
condition3 = 2 * X - 3 * Y <= 6  # для 2x - 3y <= 6
condition4 = X >= 0  # для x >= 0
condition5 = Y >= 0  # для y >= 0

# Спільна область: всі умови мають бути виконані
valid_area = condition1 & condition2 & condition3 & condition4 & condition5

# Функція для знаходження перетину з осями
def find_intersections_with_axes():
    intercepts = []
    
    # Перетини для x - y = -3
    intercepts.append(np.array([[-3, 0], [0, 3]]))  # x = -3, y = 3
    
    # Перетини для 6x + 7y = 42
    x_intercept = 42 / 6
    y_intercept = 42 / 7
    intercepts.append(np.array([[x_intercept, 0], [0, y_intercept]]))
    
    # Перетини для 2x - 3y = 6
    x_intercept = 6 / 2
    y_intercept = -2
    intercepts.append(np.array([[x_intercept, 0], [0, y_intercept]]))
    
    return intercepts

# Функція для чисельного розв'язання системи рівнянь (перетини між лініями)
def solve_system_of_equations(coeffs, constants):
    # Розв'язуємо систему лінійних рівнянь Ax = b
    return np.linalg.solve(coeffs, constants)

# Функція для знаходження точок перетину між прямими
def find_intersections_between_lines():
    intersections = []
    
    # Розв'язок системи рівнянь для x - y = -3 і 6x + 7y = 42
    A1 = np.array([[1, -1], [6, 7]])
    b1 = np.array([-3, 42])
    intersection1 = solve_system_of_equations(A1, b1)
    
    # Розв'язок системи рівнянь для 6x + 7y = 42 і 2x - 3y = 6
    A2 = np.array([[6, 7], [2, -3]])
    b2 = np.array([42, 6])
    intersection2 = solve_system_of_equations(A2, b2)
    
    # Розв'язок системи рівнянь для 2x - 3y = 6 і x - y = -3
    A3 = np.array([[2, -3], [1, -1]])
    b3 = np.array([6, -3])
    intersection3 = solve_system_of_equations(A3, b3)
    
    intersections.append(intersection1)
    intersections.append(intersection2)
    intersections.append(intersection3)
    
    return intersections

# Цільова функція F = 3x - 2y
def target_function(x, y):
    return 3 * x - 2 * y

# Перевірка, чи точка задовольняє всі умови
def is_valid_vertex(x, y):
    # Можна дозволити точкам бути на межах
    epsilon = 1e-6  # Дозволений допуск для перевірки на межах
    return (x - y >= -3 - epsilon) and (6 * x + 7 * y <= 42 + epsilon) and (2 * x - 3 * y <= 6 + epsilon) and (x >= 0) and (y >= 0)

# Знаходимо перетини
intersections_axes = find_intersections_with_axes()
intersections_lines = find_intersections_between_lines()

# Додаємо точку початку координат (0, 0) як можливу вершину
valid_vertices = [(0, 0)] + [tuple(point) for group in intersections_axes for point in group] + intersections_lines

# Фільтруємо тільки валідні точки
valid_vertices = [point for point in valid_vertices if is_valid_vertex(point[0], point[1])]

# Обчислюємо значення цільової функції на кожній точці
values_at_vertices = [target_function(x, y) for x, y in valid_vertices]

# Визначимо максимальне та мінімальне значення цільової функції
max_value = max(values_at_vertices)
min_value = min(values_at_vertices)

# Знаходимо координати точок для максимального і мінімального значення цільової функції
max_index = values_at_vertices.index(max_value)
min_index = values_at_vertices.index(min_value)
max_point = valid_vertices[max_index]
min_point = valid_vertices[min_index]

# Показати результати
print("Точки перетину з осями:")
for intercept in intersections_axes:
    print(f"Перетини: ({intercept[0][0]}, {intercept[1][1]}) і ({intercept[1][0]}, {intercept[0][1]})")

print("\nТочки перетину між лініями:")
for point in intersections_lines:
    print(f"Перетин: ({point[0]:.2f}; {point[1]:.2f})")

print("\nВершини області:")
for point in valid_vertices:
    print(f"Точка: ({point[0]:.2f}; {point[1]:.2f})")

print("\nЗначення цільової функції на цих точках:")
for (point, value) in zip(valid_vertices, values_at_vertices):
    print(f"Точка ({point[0]:.2f}; {point[1]:.2f}) має значення F = {value:.2f}")

print(f"\nМаксимум цільової функції: {max_value:.2f} при точці ({max_point[0]:.2f}; {max_point[1]:.2f})")
print(f"Мінімум цільової функції: {min_value:.2f} при точці ({min_point[0]:.2f}; {min_point[1]:.2f})")

# Побудова графіку
plt.figure(figsize=(6,6))

# Заштриховування спільної області
plt.imshow(valid_area, extent=[0, 10, 0, 10], origin="lower", cmap="Greys", alpha=0.5)

# Побудова ліній для кожної умови
plt.plot(x, x + 3, color='red', linewidth=2, linestyle='--', label=r'$x - y = -3$')
plt.plot(x, (42 - 6 * x) / 7, color='blue', linewidth=2, linestyle='--', label=r'$6x + 7y = 42$')
plt.plot(x, (2 * x - 6) / 3, color='green', linewidth=2, linestyle='--', label=r'$2x - 3y = 6$')

# Фільтруємо тільки валідні точки (точки, які задовольняють всі умови)
valid_intersections = [point for point in intersections_lines if is_valid_vertex(point[0], point[1])]

# Додаємо їх до списку валідних точок
valid_vertices.extend(valid_intersections)

# Відображення точок перетину з підписами
for i, point in enumerate(valid_vertices):
    plt.scatter(point[0], point[1], color='black')
    plt.text(point[0] + 0.1, point[1] + 0.1, f"({point[0]:.2f}; {point[1]:.2f})", fontsize=9)

# Налаштування графіку
plt.xlim(0, 10)
plt.ylim(0, 10)
plt.axhline(0, color='black',linewidth=0.7)
plt.axvline(0, color='black',linewidth=0.7)
plt.grid(True)
plt.legend()
plt.title(r'Спільна область, де виконуються всі умови')

# Показати графік
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations

def find_intersections(constraints):
    intersections = []
    for (i, j) in combinations(range(len(constraints)), 2):
        A = constraints[[i, j], :2]
        b = constraints[[i, j], 2]
        try:
            intersection = np.linalg.solve(A, b)
            intersections.append((i, j, intersection))
        except np.linalg.LinAlgError:
            pass
    return intersections

def find_axis_intersections(intersections):
    intersections_with_axes = {}
    intersections_between_lines = []
    for i, j, point in intersections:
        point = np.round(point, 6)
        if point[0] == 0 or point[1] == 0:
            axis = "Ox" if point[1] == 0 else "Oy"
            if (i, axis) not in intersections_with_axes:
                intersections_with_axes[(i, axis)] = []
            intersections_with_axes[(i, axis)].append(point)
        else:
            intersections_between_lines.append(point)
    return intersections_with_axes, intersections_between_lines

def filter_valid_vertices(intersections, constraints):
    valid_vertices = []
    for _, _, vertex in intersections:
        if np.all(np.dot(constraints[:, :2], vertex) <= constraints[:, 2] + 1e-6):
            valid_vertices.append(vertex)
    return np.array(valid_vertices)

def evaluate_objective_function(valid_vertices):
    values = []
    for vertex in valid_vertices:
        numerator = -4 * vertex[0] + 2 * vertex[1]
        denominator = -2 * vertex[0] + 3 * vertex[1]
        if abs(denominator) < 1e-8:
            values.append(np.nan)
        else:
            values.append(numerator / denominator)
    values = np.array(values)
    finite_values = np.isfinite(values)

    if np.any(finite_values):
        max_index = np.nanargmax(values)
        min_index = np.nanargmin(values)
    else:
        max_index = min_index = -1

    return values, max_index, min_index

def plot_feasible_region(constraints, valid_vertices):
    x = np.linspace(-5, 10, 400)
    y = np.linspace(-5, 10, 400)
    X, Y = np.meshgrid(x, y)
    plt.figure(figsize=(7, 7))

    colors = ['red', 'blue', 'green', 'black', 'black']
    for i, (a, b, c) in enumerate(constraints):
        label = f'{a}x + {b}y ≤ {c}' if a and b else f'{b}y ≤ {c}' if a == 0 else f'{a}x ≤ {c}'
        if b != 0:
            plt.plot(x, (c - a * x) / b, linestyle='-', label=label, color=colors[i % len(colors)])
        else:
            plt.axvline(c / a, linestyle='-', label=label, color=colors[i % len(colors)])

    valid_area = np.ones_like(X, dtype=bool)
    for (a, b, c) in constraints:
        valid_area &= (a * X + b * Y <= c)
    plt.contourf(X, Y, valid_area, levels=[0.5, 1], colors=['yellow'], alpha=0.5)

    for vertex in valid_vertices:
        plt.scatter(vertex[0], vertex[1], color='black', zorder=3)
        plt.text(vertex[0] + 0.2, vertex[1] + 0.2, f'({vertex[0]:.2f}, {vertex[1]:.2f})', fontsize=9, color='black')

    plt.axhline(0, color='black', linewidth=3)
    plt.axvline(0, color='black', linewidth=3)
    plt.xlim(-5, 10)
    plt.ylim(-5, 10)
    plt.grid(True)
    plt.legend()
    plt.title("Множина допустимих розв'язків")
    plt.show()

def main():
    constraints = np.array([
        [2, -1, 8],
        [1, 1, 5],
        [3, -2, -3],
        [-1, 0, 0],
        [0, -1, 0]
    ])

    intersections = find_intersections(constraints)
    intersections_with_axes, intersections_between_lines = find_axis_intersections(intersections)
    valid_vertices = filter_valid_vertices(intersections, constraints)

    if valid_vertices.shape[0] == 0:
        print("Помилка: Немає допустимих точок, які задовольняють всі обмеження!")
        return

    values_at_vertices, max_index, min_index = evaluate_objective_function(valid_vertices)

    print("Точки перетину з осями:")
    for (i, axis), points in intersections_with_axes.items():
        formatted_points = [f"({0.0 if round(p[0], 2) == -0.0 else round(p[0], 2)}; {0.0 if round(p[1], 2) == -0.0 else round(p[1], 2)})" for p in points]
        print(f"Пряма {i + 1} та {axis}: {', '.join(formatted_points)}")

    print("\nТочки перетину між лініями:")
    for vertex in intersections_between_lines:
        print(f"Перетин: ({vertex[0]:.2f}; {vertex[1]:.2f})")

    print("\nВершини області:")
    for vertex in valid_vertices:
        print(f"Точка: ({vertex[0]:.2f}; {vertex[1]:.2f})")

    print("\nЗначення цільової функції на цих точках:")
    for vertex, value in zip(valid_vertices, values_at_vertices):
        if np.isfinite(value):
            print(f"Точка ({vertex[0]:.2f}; {vertex[1]:.2f}) має значення F = {value:.2f}")
        else:
            print(f"Точка ({vertex[0]:.2f}; {vertex[1]:.2f}) має значення F = невизначене (ділення на 0)")

    if max_index != -1:
        max_value = values_at_vertices[max_index]
        min_value = values_at_vertices[min_index]

        max_points = [vertex for vertex, val in zip(valid_vertices, values_at_vertices)
                      if np.isclose(val, max_value, atol=1e-6)]
        min_points = [vertex for vertex, val in zip(valid_vertices, values_at_vertices)
                      if np.isclose(val, min_value, atol=1e-6)]

        max_coords = ' та при точці '.join([f'({v[0]:.2f}; {v[1]:.2f})' for v in max_points])
        min_coords = ' та при точці '.join([f'({v[0]:.2f}; {v[1]:.2f})' for v in min_points])

        print(f"\nМаксимум цільової функції: {max_value:.2f} при точці {max_coords}")
        print(f"Мінімум цільової функції: {min_value:.2f} при точці {min_coords}\n")
    else:
        print("\nУсі значення цільової функції є невизначеними (усюди ділення на 0).")

    plot_feasible_region(constraints, valid_vertices)

if __name__ == "__main__":
    main()

print("=" * 40)
print("Група: ОІ-23")
print("Виконав: Гончаров Вадим Олегович")
print("Номер ЛР: №6")
print("Варіант: №75")
print("=" * 40)
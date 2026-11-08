import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations

def clean(value):
    return 0.0 if np.isclose(value, 0) else round(value, 2)

def find_valid_vertices(constraints):
    intersections = []
    for i, j in combinations(range(len(constraints)), 2):
        A, b = constraints[[i, j], :2], constraints[[i, j], 2]
        try:
            point = np.linalg.solve(A, b)
            if np.all(constraints[:, :2] @ point <= constraints[:, 2] + 1e-6):
                intersections.append(point)
        except np.linalg.LinAlgError:
            continue
    return np.array(intersections)

def evaluate_function(points, c):
    values = points @ c
    return values, np.argmax(values), np.argmin(values)

def plot_region(constraints, vertices):
    x = np.linspace(-1, 10, 400)
    y = np.linspace(-1, 10, 400)
    X, Y = np.meshgrid(x, y)
    area = np.ones_like(X, dtype=bool)

    plt.figure(figsize=(7, 7))
    for a, b, c in constraints:
        if b != 0:
            plt.plot(x, (c - a * x) / b, linestyle='--', linewidth=2, color='black')
        else:
            plt.axvline(x=c / a, linestyle='--', color='black', linewidth=2)

    for (a, b, c) in constraints:
        area &= a * X + b * Y <= c
    plt.contourf(X, Y, area, levels=[0.5, 1], colors='gray', alpha=0.5)

    for v in vertices:
        plt.scatter(*v, color='black', zorder=3)
        x_val, y_val = clean(v[0]), clean(v[1])
        plt.text(x_val + 0.2, y_val + 0.2, f'({x_val:.2f}, {y_val:.2f})', fontsize=9)

    plt.axhline(0, color='black', linewidth=3)
    plt.axvline(0, color='black', linewidth=3)
    plt.xlim(-1, 10); plt.ylim(-1, 10)
    plt.grid(True); plt.title("Багатокутник розв'язків")
    plt.show()

def main():
    constraints = np.array([
        [-1, -2, -5],
        [2, 4, 16],
        [-6, -1, -6],
        [1, 3, 9],
        [-1, 0, 0],
        [0, -1, 0]
    ])
    c = np.array([3, 1])
    vertices = find_valid_vertices(constraints)

    if len(vertices) == 0:
        print("Немає допустимих точок!")
        return

    values, max_i, min_i = evaluate_function(vertices, c)

    print("Вершини області:")
    for v in vertices:
        print(f"({clean(v[0]):.2f}; {clean(v[1]):.2f})")

    print("\nЗначення F у вершинах:")
    for v, f in zip(vertices, values):
        print(f"F = {round(f, 2):.2f} у точці ({clean(v[0]):.2f}; {clean(v[1]):.2f})")

    vmax, vmin = round(values[max_i], 2), round(values[min_i], 2)
    vmax_point, vmin_point = vertices[max_i], vertices[min_i]

    print(f"\nmax(F) = {vmax:.2f} при ({clean(vmax_point[0]):.2f}; {clean(vmax_point[1]):.2f})")
    print(f"min(F) = {vmin:.2f} при ({clean(vmin_point[0]):.2f}; {clean(vmin_point[1]):.2f})\n")

    plot_region(constraints, vertices)

if __name__ == "__main__":
    main()
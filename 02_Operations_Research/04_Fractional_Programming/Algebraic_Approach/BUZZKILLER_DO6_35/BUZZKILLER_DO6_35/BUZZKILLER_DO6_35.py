import numpy as np
from itertools import combinations
from fractions import Fraction
from decimal import Decimal, getcontext
from tabulate import tabulate
import matplotlib.pyplot as plt

getcontext().prec = 20

def analyze_fractional_problem():
    constraints = np.array([
        [1, 1, 6],
        [3, -2, -3],
        [-1, 0, 0],
        [0, -1, 0]
    ], dtype=float)

    vertices = calculate_intersection_points(constraints)
    feasible = filter_feasible_vertices(vertices, constraints)
    if not len(feasible):
        print("Допустимих розв’язків не знайдено.")
        return

    f_values, idx_max, idx_min = evaluate_fractional_function(feasible)
    display_result_table(feasible, f_values, idx_max, idx_min)
    plot_feasible_area_and_points(constraints, feasible)

def calculate_intersection_points(ineqs):
    result = []
    for i, j in combinations(range(len(ineqs)), 2):
        A, b = ineqs[[i, j], :2], ineqs[[i, j], 2]
        try:
            solution = np.linalg.solve(A, b)
            result.append((i, j, solution))
        except np.linalg.LinAlgError:
            continue
    return result

def filter_feasible_vertices(points, constraints):
    return np.array([pt for _, _, pt in points if np.all(constraints[:, :2] @ pt <= constraints[:, 2] + 1e-8)])

def evaluate_fractional_function(vertices):
    results = []
    for v in vertices:
        x1, x2 = Decimal(str(v[0])), Decimal(str(v[1]))
        num = -2 * x1 + x2
        den = 2 * x1 - 3 * x2
        if abs(den) < Decimal("1e-8"):
            results.append(None)
        else:
            results.append(Fraction(num / den).limit_denominator())

    idx_max = idx_min = -1
    finite = [(i, val) for i, val in enumerate(results) if val is not None]
    if finite:
        max_val = max(v for _, v in finite)
        min_val = min(v for _, v in finite)
        idx_max = next(i for i, v in finite if v == max_val)
        idx_min = next(i for i, v in finite if v == min_val)
    return results, idx_max, idx_min

def display_result_table(vertices, values, idx_max, idx_min):
    to_frac = lambda val: Fraction(Decimal(str(val))).limit_denominator()
    rows = [[f"({to_frac(p[0])}; {to_frac(p[1])})", str(v) if v is not None else "невизначене"]
            for p, v in zip(vertices, values)]
    print(tabulate(rows, headers=["Допустима вершина", "F(x₁, x₂)"], tablefmt="simple"))

    if idx_max == -1:
        print("\nF(x₁, x₂) невизначена для всіх вершин.")
        return

    max_val = values[idx_max]
    min_val = values[idx_min]
    max_pts = [f"({to_frac(p[0])}; {to_frac(p[1])})" for p, v in zip(vertices, values) if v == max_val]
    min_pts = [f"({to_frac(p[0])}; {to_frac(p[1])})" for p, v in zip(vertices, values) if v == min_val]

    print("\n" + tabulate([
        ["Максимум", str(max_val), ", ".join(max_pts)],
        ["Мінімум", str(min_val), ", ".join(min_pts)],
    ], headers=["Екстремум", "Значення", "Досягається в точках"], tablefmt="simple"))

def plot_feasible_area_and_points(bounds, vertices):
    x = np.linspace(-1, 10, 800)
    y = np.linspace(-1, 10, 800)
    X, Y = np.meshgrid(x, y)
    mask = np.ones_like(X, dtype=bool)

    plt.figure(figsize=(10, 7))
    for a, b, c in bounds:
        mask &= a * X + b * Y <= c
        if b != 0:
            plt.plot(x, (c - a * x) / b, '--', lw=1.5)
        else:
            plt.axvline(x=c / a, linestyle='--', lw=1.5)

    plt.contourf(X, Y, mask, levels=[0.5, 1], colors=['blue'], alpha=0.4)

    for pt in vertices:
        x1f = Fraction(Decimal(str(pt[0]))).limit_denominator()
        x2f = Fraction(Decimal(str(pt[1]))).limit_denominator()
        plt.scatter(pt[0], pt[1], color='red')
        plt.text(pt[0] + 0.1, pt[1] + 0.1, f"({x1f}; {x2f})", fontsize=9)

    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.xlim(-1, 10)
    plt.ylim(-1, 10)
    plt.xlabel("x₁")
    plt.ylabel("x₂")
    plt.title("Область допустимих розв’язків")
    plt.grid(True, linestyle=':')
    plt.tight_layout()
    plt.show()

analyze_fractional_problem()

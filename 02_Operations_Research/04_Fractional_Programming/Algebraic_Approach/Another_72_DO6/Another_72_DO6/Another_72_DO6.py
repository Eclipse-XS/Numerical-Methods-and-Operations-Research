import numpy as np
from itertools import combinations
from fractions import Fraction
from decimal import Decimal, getcontext
from tabulate import tabulate
import matplotlib.pyplot as plt

getcontext().prec = 20

def analyze_constraints(ineqs):
    points = []
    for i, j in combinations(range(len(ineqs)), 2):
        A, b = ineqs[[i, j], :2], ineqs[[i, j], 2]
        try:
            sol = np.linalg.solve(A, b)
            points.append((i, j, sol))
        except np.linalg.LinAlgError:
            continue
    return points

def filter_feasible_points(points, ineqs):
    return np.array([pt for _, _, pt in points if np.all(ineqs[:, :2] @ pt <= ineqs[:, 2] + 1e-8)])

def evaluate_fractional_function(vertices):
    values = []
    for v in vertices:
        x1, x2 = Decimal(str(v[0])), Decimal(str(v[1]))
        num = Decimal(2) * x1 + Decimal(2) * x2
        den = Decimal(-2) * x1 + Decimal(3) * x2
        if abs(den) < Decimal("1e-8"):
            values.append(None)
        else:
            values.append(Fraction(num / den).limit_denominator())
    idx_max = idx_min = -1
    finite_vals = [(i, v) for i, v in enumerate(values) if v is not None]
    if finite_vals:
        max_val = max(v for _, v in finite_vals)
        min_val = min(v for _, v in finite_vals)
        idx_max = next(i for i, v in finite_vals if v == max_val)
        idx_min = next(i for i, v in finite_vals if v == min_val)
    return values, idx_max, idx_min

def display_table_results(vertices, values, idx_max, idx_min):
    to_frac = lambda val: Fraction(Decimal(str(val))).limit_denominator()
    table = [[f"({to_frac(x[0])}; {to_frac(x[1])})", str(v) if v is not None else "невизначене"]
             for x, v in zip(vertices, values)]
    print(tabulate(table, headers=["Допустима вершина", "F(x1, x2)"], tablefmt="github"))

    if idx_max == -1:
        print("\nF(x1, x2) невизначена для всіх точок.")
        return

    max_val, min_val = values[idx_max], values[idx_min]
    max_pts = [f"({to_frac(x[0])}; {to_frac(x[1])})"
               for x, v in zip(vertices, values) if v == max_val]
    min_pts = [f"({to_frac(x[0])}; {to_frac(x[1])})"
               for x, v in zip(vertices, values) if v == min_val]

    print("\n" + tabulate([
        ["Максимум", str(max_val), ", ".join(max_pts)],
        ["Мінімум", str(min_val), ", ".join(min_pts)],
    ], headers=["Екстремум", "Значення", "У точках"], tablefmt="github"))

def plot_feasible_area(bounds, vertices):
    x = np.linspace(-2, 10, 800)
    y = np.linspace(-2, 10, 800)
    X, Y = np.meshgrid(x, y)
    mask = np.ones_like(X, dtype=bool)

    plt.figure(figsize=(9, 7))
    for a, b, c in bounds:
        mask &= a * X + b * Y <= c
        if b != 0:
            plt.plot(x, (c - a * x) / b, lw=1.5, ls='--')
        else:
            plt.axvline(x=c / a, lw=1.5, ls='--')

    plt.contourf(X, Y, mask, levels=[0.5, 1], colors=['#cce5ff'], alpha=0.4)

    for pt in vertices:
        x1f = Fraction(Decimal(str(pt[0]))).limit_denominator()
        x2f = Fraction(Decimal(str(pt[1]))).limit_denominator()
        plt.scatter(pt[0], pt[1], color='red')
        plt.text(pt[0] + 0.1, pt[1] + 0.1, f"({x1f}; {x2f})", fontsize=8)

    plt.axhline(0, color='black', lw=1)
    plt.axvline(0, color='black', lw=1)
    plt.xlim(-2, 10)
    plt.ylim(-2, 10)
    plt.xlabel("x₁")
    plt.ylabel("x₂")
    plt.title("Область допустимих розв'язків")
    plt.grid(True, linestyle=':')
    plt.tight_layout()
    plt.show()

def run_analysis():
    constraints = np.array([
        [-1, -1, -1],    # x1 + x2 ≥ 1
        [-5,  1,  0],    # -5x1 + x2 ≤ 0
        [ 1, -5,  0],    # -x1 + 5x2 ≥ 0
        [ 1,  1,  6],    # x1 + x2 ≤ 6
        [-1,  0,  0],    # x1 ≥ 0
        [ 0, -1,  0],    # x2 ≥ 0
    ], dtype=float)

    points = analyze_constraints(constraints)
    feasible = filter_feasible_points(points, constraints)
    if not len(feasible):
        print("Немає допустимих точок!")
        return

    values, idx_max, idx_min = evaluate_fractional_function(feasible)
    display_table_results(feasible, values, idx_max, idx_min)
    plot_feasible_area(constraints, feasible)

if __name__ == "__main__":
    run_analysis()
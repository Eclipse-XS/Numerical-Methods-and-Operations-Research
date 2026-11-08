import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations

def get_line_intersections(bounds):
    pts = []
    for i, j in combinations(range(len(bounds)), 2):
        A, b = bounds[[i, j], :2], bounds[[i, j], 2]
        try:
            pts.append((i, j, np.linalg.solve(A, b)))
        except np.linalg.LinAlgError:
            continue
    return pts

def separate_points(pts):
    on_axes, inner = {}, []
    for i, j, pt in pts:
        pt = np.round(pt, 6)
        if pt[0] == 0 or pt[1] == 0:
            axis = 'Ox' if pt[1] == 0 else 'Oy'
            on_axes.setdefault((i, axis), []).append(pt)
        else:
            inner.append(pt)
    return on_axes, inner

def valid_vertices(pts, bounds):
    return np.array([pt for _, _, pt in pts if np.all(bounds[:, :2] @ pt <= bounds[:, 2] + 1e-6)])

def fractional_eval(vertices):
    vals = []
    for v in vertices:
        num = 2 * v[0] + 2 * v[1]
        den = -2 * v[0] + 3 * v[1]
        if abs(den) < 1e-8:
            vals.append(np.nan)
        else:
            vals.append(num / den)
    vals = np.array(vals)
    return vals, np.nanargmax(vals) if np.any(np.isfinite(vals)) else -1, np.nanargmin(vals) if np.any(np.isfinite(vals)) else -1

def summary_output(axes, between, vertices, values, idx_max, idx_min):
    print("\n[Перетини з осями]")
    for (i, axis), pts in axes.items():
        print(f"Пряма {i+1} ∩ {axis}:", ', '.join(f"({p[0]:.2f}; {p[1]:.2f})" for p in pts))
    print("\n[Перетини між прямими]")
    for pt in between:
        print(f"({pt[0]:.2f}; {pt[1]:.2f})")
    print("\n[Вершини області допуску]")
    for pt in vertices:
        print(f"({pt[0]:.2f}; {pt[1]:.2f})")
    print("\n[Значення цільової функції]")
    for pt, v in zip(vertices, values):
        res = f"{v:.2f}" if np.isfinite(v) else "невизначене"
        print(f"F({pt[0]:.2f}; {pt[1]:.2f}) = {res}")
    if idx_max == -1:
        return print("\nF невизначена всюди")
    def ext_pts(val): return [p for p, v in zip(vertices, values) if np.isclose(v, val)]
    vmax, vmin = values[idx_max], values[idx_min]
    max_pts, min_pts = ext_pts(vmax), ext_pts(vmin)
    print(f"\nМаксимум F = {vmax:.2f} у точках {', '.join(f'({p[0]:.2f}; {p[1]:.2f})' for p in max_pts)}")
    print(f"Мінімум F = {vmin:.2f} у точках {', '.join(f'({p[0]:.2f}; {p[1]:.2f})' for p in min_pts)}")

def draw_region(bounds, vertices):
    x = np.linspace(-1, 10, 400)
    y = np.linspace(-1, 10, 400)
    X, Y = np.meshgrid(x, y)
    plt.figure(figsize=(8, 6))
    mask = np.ones_like(X, dtype=bool)
    for a, b, c in bounds:
        mask &= a * X + b * Y <= c
        if b != 0:
            plt.plot(x, (c - a * x) / b, lw=2, ls='--')
        else:
            plt.axvline(x=c / a, lw=2, ls='--')
    plt.contourf(X, Y, mask, levels=[0.5, 1], colors=['#d1ecf1'], alpha=0.6)
    for pt in vertices:
        plt.scatter(*pt, color='black')
        plt.text(pt[0] + 0.1, pt[1] + 0.1, f"({pt[0]:.1f}; {pt[1]:.1f})", fontsize=8)
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.xlim(-1, 10)
    plt.ylim(-1, 10)
    plt.grid(True, linestyle=':')
    plt.title("Область допустимих розв'язків")
    plt.tight_layout()
    plt.show()

def main():
    bounds = np.array([
        [-1, -1, -1],     # x1 + x2 ≥ 1     → -x1 - x2 ≤ -1
        [-5,  1,  0],     # -5x1 + x2 ≤ 0
        [ 1, -5,  0],     # -x1 + 5x2 ≥ 0   → x1 - 5x2 ≤ 0
        [ 1,  1,  6],     # x1 + x2 ≤ 6
        [-1,  0,  0],     # x1 ≥ 0          → -x1 ≤ 0
        [ 0, -1,  0],     # x2 ≥ 0          → -x2 ≤ 0
    ])
    pts = get_line_intersections(bounds)
    axes, inner = separate_points(pts)
    feasible = valid_vertices(pts, bounds)
    if not len(feasible):
        return print("Немає допустимих точок!")
    vals, idx_max, idx_min = fractional_eval(feasible)
    summary_output(axes, inner, feasible, vals, idx_max, idx_min)
    draw_region(bounds, feasible)

if __name__ == "__main__":
    main()

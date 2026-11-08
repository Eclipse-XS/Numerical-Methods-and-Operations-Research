import sympy as sp
from sympy import Matrix, simplify
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations

# ---------------------- ГРАФІЧНИЙ МЕТОД ----------------------

def round_val(val):
    return 0.0 if np.isclose(val, 0) else round(val, 2)

def get_feasible_vertices(ineqs):
    pts = []
    for i, j in combinations(range(len(ineqs)), 2):
        A = ineqs[[i, j], :2]
        b = ineqs[[i, j], 2]
        try:
            sol = np.linalg.solve(A, b)
            if np.all(ineqs[:, :2] @ sol <= ineqs[:, 2] + 1e-6):
                pts.append(tuple(np.round(sol, 6)))
        except np.linalg.LinAlgError:
            continue
    return np.unique(np.array(pts), axis=0)

def fractional_value(points):
    results = []
    for x1, x2 in points:
        num = 2 * x1 + 3 * x2
        den = 2 * x1 + x2
        results.append(np.inf if np.isclose(den, 0) else num / den)
    results = np.array(results)
    return results, np.where(np.isclose(results, np.max(results)))[0], np.where(np.isclose(results, np.min(results)))[0]

def draw_region(ineqs, verts):
    x = np.linspace(-1, 6, 800)
    y = np.linspace(-1, 6, 800)
    X, Y = np.meshgrid(x, y)
    mask = np.ones_like(X, dtype=bool)

    plt.figure(figsize=(8, 8))
    for a, b, c in ineqs:
        if b != 0:
            plt.plot(x, (c - a * x) / b, '--', color='black', linewidth=1.5)
        else:
            plt.axvline(x=c / a, linestyle='--', color='black', linewidth=1.5)

    for a, b, c in ineqs:
        mask &= a * X + b * Y <= c
    plt.contourf(X, Y, mask, levels=[0.5, 1], colors='skyblue', alpha=0.6)

    for x1, x2 in verts:
        plt.scatter(x1, x2, color='black', zorder=3)
        plt.text(x1 + 0.1, x2 + 0.1, f'({round_val(x1):.2f}; {round_val(x2):.2f})', fontsize=9)

    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.xlim(-1, 6)
    plt.ylim(-1, 6)
    plt.grid(True)
    plt.title("Область допустимих розв’язків")
    plt.show()

def solve_graphically():
    constraints = np.array([
        [-3, 2, 6],
        [-1, -2, -3],
        [1, 0, 3],
        [-1, 0, 0],
        [0, 1, 5],
        [0, -1, 0]
    ])

    vertices = get_feasible_vertices(constraints)
    if len(vertices) == 0:
        print("Немає допустимих точок.")
        return

    print("Вершини області допустимих розв’язків:")
    for x1, x2 in vertices:
        print(f"({round_val(x1):.2f}; {round_val(x2):.2f})")

    values, max_idx, min_idx = fractional_value(vertices)

    print("\nМаксимум F(x):")
    for i in max_idx:
        x1, x2 = round_val(vertices[i][0]), round_val(vertices[i][1])
        print(f"F = {round(values[i], 4):.4f} у точці ({x1:.2f}; {x2:.2f})")

    print("\nМінімум F(x):")
    for i in min_idx:
        x1, x2 = round_val(vertices[i][0]), round_val(vertices[i][1])
        print(f"F = {round(values[i], 4):.4f} у точці ({x1:.2f}; {x2:.2f})")

    draw_region(constraints, vertices)

# ---------------------- СИМПЛЕКС-МЕТОД ----------------------

def linearize(num, den, A, b, signs, goal):
    num = list(map(sp.Rational, num))
    den = list(map(sp.Rational, den))
    A = [[sp.Rational(a) for a in row] for row in A]
    b = list(map(sp.Rational, b))
    n = len(num)
    obj = num + [0]
    A_lin, b_lin, signs_lin = [], [], []
    for i in range(len(A)):
        if b[i] < 0:
            A[i] = [-a for a in A[i]]
            b[i] = -b[i]
            signs[i] = ">=" if signs[i] == "<=" else "<=" if signs[i] == ">=" else signs[i]
        row = A[i] + [-b[i]]
        if signs[i] == ">=":
            row = [-a for a in row]
            signs_lin.append("<=")
            b_lin.append(0)
        else:
            signs_lin.append(signs[i])
            b_lin.append(0)
        A_lin.append(row)
    A_lin.append(den + [0])
    signs_lin.append("=")
    b_lin.append(1)
    return A_lin, obj, [f"y{i+1}" for i in range(n)] + ["y0"], signs_lin, b_lin

def build_tableau(obj, A, b, signs, goal):
    A, b = Matrix(A), Matrix(b).reshape(len(b), 1)
    T = A
    names = [f"y{i+1}" for i in range(A.shape[1] - 1)] + ["y0"]
    base, art = [], []
    s_id = a_id = 1
    for i, s in enumerate(signs):
        if s == "<=":
            T = T.row_join(Matrix.eye(len(signs)).col(i))
            v = f"s{s_id}"; names.append(v); base.append(v); s_id += 1
        elif s == ">=":
            T = T.row_join(-Matrix.eye(len(signs)).col(i))
            T = T.row_join(Matrix.eye(len(signs)).col(i))
            names += [f"s{s_id}", f"w{a_id}"]
            base.append(f"w{a_id}"); art.append((i, f"w{a_id}")); s_id += 1; a_id += 1
        elif s == "=":
            T = T.row_join(Matrix.eye(len(signs)).col(i))
            v = f"w{a_id}"; names.append(v); base.append(v); art.append((i, v)); a_id += 1
    T = T.row_join(b)
    c = Matrix([sp.Rational(obj[i]) if i < len(obj) else 0 for i in range(T.cols)])
    mult = 1 if goal == 'min' else -1
    for idx, av in art:
        row = T[idx, :]
        for j in range(len(c) - 1):
            if names[j] != av: c[j] += mult * M * (-row[j])
        c[names.index(av)] = 0; c[-1] += mult * M * row[-1]
    for i in range(len(c) - 1): c[i] *= -1
    return T, c, names, base

def simplex_solver(T, c, names, base, goal):
    while True:
        row = [float(v.subs(M, 1000)) for v in c[:-1]]
        col = (row.index(min(row)) if min(row) < 0 else None) if goal == 'max' else (
            max((i for i, v in enumerate(row) if v > 0), key=lambda i: row[i], default=None))
        if col is None: break
        ratios = [(T[i, -1] / T[i, col] if T[i, col] > 0 else float('inf')) for i in range(T.shape[0])]
        if all(r == float('inf') for r in ratios): return None, None, None
        r = ratios.index(min(ratios, key=float))
        e = T[r, col]
        T.row_op(r, lambda x, _: simplify(x / e))
        for i in range(T.shape[0]):
            if i != r:
                f = T[i, col]
                T.row_op(i, lambda x, j: simplify(x - f * T[r, j]))
        f = c[col]
        c = Matrix([simplify(c[j] - f * T[r, j]) for j in range(c.shape[0])])
        base[r] = names[col]
    sol = [0] * c.shape[0]
    for i in range(T.shape[0]):
        sol[names.index(base[i])] = T[i, -1]
    return sol, c[-1], base

def solve_by_simplex():
    for goal in ["max", "min"]:
        print(f"\n{'='*10} РОЗВ'ЯЗОК ДЛЯ {goal.upper()} {'='*10}")
        num = [2, 3]
        den = [2, 1]
        A = [[3, -2], [1, 2], [1, 0], [0, 1]]
        b = [-6, 3, 3, 5]
        signs = [">=", ">=", "<=", "<="]
        A_lin, c_lin, vars_, s_new, b_lin = linearize(num, den, A, b, signs, goal)
        T, c_row, names, base = build_tableau(c_lin, A_lin, b_lin, s_new, goal)
        sol, val, _ = simplex_solver(T.copy(), c_row.copy(), names, base, goal)
        if sol:
            y0 = sol[names.index("y0")]
            if y0 == 0:
                print("y0 = 0 → зворотна заміна неможлива (ділення на нуль).")
            else:
                xs = [simplify(sol[names.index(f"y{i+1}")] / y0) for i in range(len(num))]
                for i, val in enumerate(xs):
                    print(f"  x{i+1} = {val}")
                numerator = sum(num[i] * xs[i] for i in range(len(xs)))
                denominator = sum(den[i] * xs[i] for i in range(len(xs)))
                f_val = simplify(numerator / denominator)
                print(f"  F(x) = {f_val}")
        else:
            print("Розв’язок не знайдено.")

# ---------------------- ВИКЛИК ----------------------

def main():
    method = input("Оберіть метод розв’язання (graph/simplex): ").strip().lower()
    if method == "graph":
        solve_graphically()
    elif method == "simplex":
        solve_by_simplex()
    else:
        print("Невідомий метод. Оберіть 'graph' або 'simplex'.")

M = sp.Symbol("M")

main()
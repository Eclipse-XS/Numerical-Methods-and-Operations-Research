import sympy as sp
from sympy import Matrix, simplify
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from fractions import Fraction
from tabulate import tabulate

M = sp.Symbol("M")

def round_value(x):
    return 0.0 if np.isclose(x, 0) else round(x, 2)

def find_vertices(constraints):
    vertices = []
    for i, j in combinations(range(len(constraints)), 2):
        A, b = constraints[[i, j], :2], constraints[[i, j], 2]
        try:
            point = np.linalg.solve(A, b)
            if np.all(constraints[:, :2] @ point <= constraints[:, 2] + 1e-6):
                if not np.isclose(point[0], 0) or not np.isclose(point[1], 0):
                    vertices.append(point)
        except np.linalg.LinAlgError:
            continue
    return np.array(vertices)

def compute_objective(points):
    values = []
    for x1, x2 in points:
        num = 3 * x1 + 4 * x2
        den = 2 * x1 + x2
        if np.isclose(den, 0):
            values.append(np.inf)
        else:
            values.append(Fraction(num).limit_denominator() / Fraction(den).limit_denominator())
    return np.array(values, dtype=object), np.argmax([float(v) for v in values]), np.argmin([float(v) for v in values])

def plot_feasible_region(constraints, vertices):
    x = np.linspace(-1, 10, 400)
    y = np.linspace(-1, 10, 400)
    X, Y = np.meshgrid(x, y)
    feasible = np.ones_like(X, dtype=bool)
    plt.figure(figsize=(7, 7))
    for a, b, c in constraints:
        if b != 0:
            plt.plot(x, (c - a * x) / b, '--', color='black')
        else:
            plt.axvline(x=c / a, linestyle='--', color='black')
    for a, b, c in constraints:
        feasible &= a * X + b * Y <= c
    plt.contourf(X, Y, feasible, levels=[0.5, 1], colors='gray', alpha=0.5)
    for v in vertices:
        plt.scatter(*v, color='black')
        plt.text(v[0] + 0.2, v[1] + 0.2, f'({round_value(v[0]):.2f}, {round_value(v[1]):.2f})', fontsize=9)
    plt.axhline(0, color='black'), plt.axvline(0, color='black')
    plt.xlim(-1, 10), plt.ylim(-1, 10), plt.grid(True)
    plt.title("Feasible Region")
    plt.show()

def solve_graphically():
    constraints = np.array([
        [8, -5, 40],
        [2, 5, 10],
        [-6, 5, 60],
        [2, 1, 14],
        [-1, 0, 0],
        [0, -1, 0]
    ])
    vertices = find_vertices(constraints)
    if len(vertices) == 0:
        print("No feasible points.")
        return
    values, max_i, min_i = compute_objective(vertices)
    table = [[f"({round_value(v[0]):.2f}; {round_value(v[1]):.2f})", str(f) if np.isfinite(float(f)) else "∞"]
             for v, f in zip(vertices, values)]
    print(tabulate(table, headers=["Vertex", "F(x1, x2)"], tablefmt="grid"))
    print(f"\nMax F(x) = {values[max_i]} at ({round_value(vertices[max_i][0])}; {round_value(vertices[max_i][1])})")
    print(f"Min F(x) = {values[min_i]} at ({round_value(vertices[min_i][0])}; {round_value(vertices[min_i][1])})")
    plot_feasible_region(constraints, vertices)

def linearize(num, den, A, b, signs, goal):
    num, den = list(map(sp.Rational, num)), list(map(sp.Rational, den))
    A = [[sp.Rational(a) for a in row] for row in A]
    b = list(map(sp.Rational, b))
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
    return A_lin, obj, [f"y{i+1}" for i in range(len(num))] + ["y0"], signs_lin, b_lin

def build_tableau(obj, A, b, signs, goal):
    A, b = Matrix(A), Matrix(b).reshape(len(b), 1)
    T, names = A, [f"y{i+1}" for i in range(A.shape[1] - 1)] + ["y0"]
    base, art, s_id, a_id = [], [], 1, 1
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
        num, den = [3, 4], [2, 1]
        A = [[8, -5], [2, 5], [-6, 5], [2, 1]]
        b = [40, 10, 60, 14]
        signs = ["<=", ">=", "<=", "<="]
        print("Дробово-лінійна задача:")
        x_vars = [sp.Symbol(f"x{i+1}") for i in range(len(num))]
        F_expr = sp.Add(*[a * x for a, x in zip(num, x_vars)]) / sp.Add(*[a * x for a, x in zip(den, x_vars)])
        print(f"F = {F_expr} → {goal}")
        print("За умов:")
        for i in range(len(A)):
            lhs = sp.Add(*[A[i][j] * x_vars[j] for j in range(len(x_vars))])
            print(lhs, signs[i], b[i])
        A_lin, c_lin, vars_, s_new, b_lin = linearize(num, den, A, b, signs, goal)
        T, c_row, names, base = build_tableau(c_lin, A_lin, b_lin, s_new, goal)
        sol, val, _ = simplex_solver(T.copy(), c_row.copy(), names, base, goal)
        print(f"Оптимальний план:")
        if sol:
            for i, v in enumerate(vars_):
                print(f" {v} = {simplify(sol[names.index(v)])}")
            y0 = sol[names.index("y0")]
            if y0 == 0:
                print("y0 = 0 → зворотна заміна неможлива (ділення на нуль).")
            else:
                print("Зворотня заміна:")
                xs = [simplify(sol[names.index(f"y{i+1}")] / y0) for i in range(len(num))]
                for i, val in enumerate(xs):
                    print(f" x{i+1} = {val}")
                numerator = sum(num[i] * xs[i] for i in range(len(xs)))
                denominator = sum(den[i] * xs[i] for i in range(len(xs)))
                print(f" F(x) = {simplify(numerator / denominator)}")
        else:
            print("Розв’язок не знайдено.")
def main():
    print("=== Варіант №91 ===")
    print("F(x1, x2) = (3*x1 + 4*x2) / (2*x1 + x2)")
    print("Умови:")
    print("8*x1 - 5*x2 ≤ 40")
    print("2*x1 + 5*x2 ≥ 10")
    print("-6*x1 + 5*x2 ≤ 60")
    print("2*x1 + x2 ≤ 14")
    print("x1 ≥ 0")
    print("x2 ≥ 0\n")

    print("Оберіть метод розв’язання:")
    print("1 — Графічний метод")
    print("2 — Лінеаризація (симплекс-метод)")
    try:
        choice = int(input("Ваш вибір (1 або 2): ").strip())
    except Exception:
        print("Некоректне введення. Очікується число 1 або 2.")
        return
    if choice == 1:
        print("=== ГРАФІЧНИЙ МЕТОД ===")
        solve_graphically()
    elif choice == 2:
        print("=== МЕТОД ЛІНЕАРИЗАЦІЇ (СИМПЛЕКС-МЕТОД) ===")
        solve_by_simplex()
    else:
        print("Невідомий метод. Оберіть 1 або 2.")


main()

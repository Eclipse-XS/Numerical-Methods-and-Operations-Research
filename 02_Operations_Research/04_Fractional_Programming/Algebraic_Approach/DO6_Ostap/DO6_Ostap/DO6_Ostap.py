
import sympy as sp
from sympy import Matrix, simplify
M = sp.Symbol("M")
G = 'min'
num = [3, 1]
den = [2, 1]
A = [[1, 2], [2, 4], [6, 1], [1, 3]]
b = [5, 16, 6, 9]
signs = [">=", "<=", ">=", "<="]

def print_function_and_constraints(obj, A, b, vars_, signs, goal):
    sym_vars = [sp.Symbol(v) for v in vars_]
    print("Лінеаризована задача:")
    print(f"G = {sp.Add(*[a*v for a,v in zip(obj, sym_vars)])} → {goal}")
    for i in range(len(A)):
        lhs = sp.Add(*[A[i][j]*sym_vars[j] for j in range(len(vars_))])
        print(f"   {lhs} {signs[i]} {b[i]}")

def build_tableau(obj, A, b, signs, goal):
    A, b = Matrix(A), Matrix(b).reshape(len(b), 1); T = A
    names = [f"y{i+1}" for i in range(A.shape[1] - 1)] + ["y0"]; base, art = [], []
    s_id = a_id = 1
    for i, s in enumerate(signs):
        if s == "<=" : T = T.row_join(Matrix.eye(len(signs)).col(i)); v = f"s{s_id}"; names.append(v); base.append(v); s_id += 1
        elif s == ">=" :
            T = T.row_join(-Matrix.eye(len(signs)).col(i)); T = T.row_join(Matrix.eye(len(signs)).col(i))
            names += [f"s{s_id}", f"w{a_id}"]; base.append(f"w{a_id}"); art.append((i, f"w{a_id}")); s_id += 1; a_id += 1
        elif s == "=":
            T = T.row_join(Matrix.eye(len(signs)).col(i)); v = f"w{a_id}"; names.append(v); base.append(v); art.append((i, v)); a_id += 1
    T = T.row_join(b); c = Matrix([sp.Rational(obj[i]) if i < len(obj) else 0 for i in range(T.cols)])
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
        r = ratios.index(min(ratios, key=float)); e = T[r, col]
        T.row_op(r, lambda x, _: simplify(x / e))
        for i in range(T.shape[0]):
            if i != r:
                f = T[i, col]; T.row_op(i, lambda x, j: simplify(x - f * T[r, j]))
        f = c[col]; c = Matrix([simplify(c[j] - f * T[r, j]) for j in range(c.shape[0])]); base[r] = names[col]
    sol = [0] * c.shape[0]
    for i in range(T.shape[0]): sol[names.index(base[i])] = T[i, -1]
    return sol, c[-1], base

def linearize(num, den, A, b, signs, goal):
    num = list(map(sp.Rational, num)); den = list(map(sp.Rational, den))
    A = [[sp.Rational(a) for a in row] for row in A]; b = list(map(sp.Rational, b))
    n = len(num); obj = num + [0]; A_lin = []; b_lin = []; signs_lin = []
    for i in range(len(A)):
        if b[i] < 0:
            A[i] = [-a for a in A[i]]; b[i] = -b[i]
            signs[i] = ">=" if signs[i] == "<=" else "<=" if signs[i] == ">=" else signs[i]
        row = A[i] + [-b[i]]
        if signs[i] == ">=": row = [-a for a in row]; signs_lin.append("<="); b_lin.append(0)
        else: signs_lin.append(signs[i]); b_lin.append(0)
        A_lin.append(row)
    A_lin.append(den + [0]); signs_lin.append("="); b_lin.append(1)
    return A_lin, obj, [f"y{i+1}" for i in range(n)] + ["y0"], signs_lin, b_lin

A_lin, c_lin, vars_, s_new, b_lin = linearize(num, den, A, b, signs, G)
print_function_and_constraints(c_lin, A_lin, b_lin, vars_, s_new, G)
T, c_row, names, base = build_tableau(c_lin, A_lin, b_lin, s_new, G)
sol, val, _ = simplex_solver(T.copy(), c_row.copy(), names, base, G)
if sol:
    print("\nЗначення змінних y:")
    for var in vars_:
        print(f"  {var} = {simplify(sol[names.index(var)])}")
    y0 = sol[names.index("y0")]
    if y0 == 0:
        print("\ny0 = 0 → зворотна заміна неможлива (ділення на нуль).")
    else:
        print("\nЗворотна заміна (xᵢ = yᵢ / y₀):")
        xs = [simplify(sol[names.index(f"y{i+1}")] / y0) for i in range(len(num))]
        for i, val in enumerate(xs):
            print(f"  x{i+1} = {val}")
        numerator = sum(num[i] * xs[i] for i in range(len(xs)))
        denominator = sum(den[i] * xs[i] for i in range(len(xs)))
        f_val = simplify(numerator / denominator)
        print("\nОптимальне значення дробово-раціональної цільової функції:")
        print(f"  F(x) = {f_val}")
else:
    print("Розв’язок не знайдено.")
import numpy as np
import sympy as sp
from sympy import Matrix, simplify

M = sp.Symbol("M")

def expr(coefs, p='x'):
    return " ".join([f"{'' if a < 0 or i == 0 else '+'}{a}{p}{i+1}" for i, a in enumerate(coefs)])

def adjust_constraints(A, b, signs, sense, flip=False):
    desired = '>=' if sense == 'min' else '<='
    A2, b2, s2 = [], [], []
    for i in range(len(b)):
        row, val, sign = A[i], b[i], signs[i]
        if (flip and val < 0) or (not flip and sign != desired):
            row, val = [-x for x in row], -val
            sign = {'>=': '<=', '<=': '>=', '=': '='}[sign]
        A2.append(row); b2.append(val); s2.append(sign)
    return A2, b2, s2

def dual(A, b, c, sense='min'):
    A, b, c = np.array(A), np.array(b), np.array(c)
    return {
        'type': 'max' if sense == 'min' else 'min',
        'A': A.T.tolist(), 'b': c.tolist(), 'c': b.tolist(),
        'signs': ['<=' if sense == 'min' else '>='] * len(A[0])
    }

def simplex_table(c, A, b, signs, goal='min'):
    A, b, vars, base, table = Matrix(A), Matrix(b).reshape(len(b), 1), [f"y{i+1}" for i in range(len(c))], [], Matrix(A)
    arts = []
    for i, s in enumerate(signs):
        if s in ['<=', '≥']:
            slack = Matrix.zeros(len(b), 1); slack[i] = 1 if s == '<=' else -1
            table = table.row_join(slack); vars.append(f"s{i+1}")
            if s == '<=': base.append(vars[-1])
        if s in ['>=', '=', '≥']:
            art = Matrix.zeros(len(b), 1); art[i] = 1
            table = table.row_join(art); vars.append(f"w{i+1}")
            base.append(vars[-1]); arts.append((i, table.cols - 1))
    table = table.row_join(b)
    z = Matrix([0]*(table.cols-1) + [0])
    for i in range(len(c)): z[i] = sp.Rational(c[i])
    for i, _ in arts:
        row, rhs = table.row(i), table[i, -1]
        for j in range(len(z)-1): z[j] += (M if goal == 'min' else -M) * (-row[j])
        z[-1] += (M if goal == 'min' else -M) * rhs
    for i in range(len(z)-1): z[i] *= -1
    for _, col in arts: z[col] = 0
    return table, z, vars, base

def solve(table, z, names, base, goal='min'):
    while True:
        numz = [float(e.subs(M, 1000)) for e in z[:-1]]
        pc = (max if goal == 'min' else min)(
            (i for i, v in enumerate(numz) if (v > 0 if goal == 'min' else v < 0)),
            default=None, key=lambda i: numz[i]
        )
        if pc is None: break
        ratios = [table[i, -1] / table[i, pc] if table[i, pc] > 0 else float('inf') for i in range(table.rows)]
        if all(r == float('inf') for r in ratios): return None, None, None
        pr = ratios.index(min(ratios)); pivot = table[pr, pc]
        table.row_op(pr, lambda x, _: simplify(x / pivot))
        for i in range(table.rows):
            if i != pr:
                factor = table[i, pc]
                table.row_op(i, lambda x, j: simplify(x - factor * table[pr, j]))
        factor = z[pc]
        z = Matrix([simplify(z[j] - factor * table[pr, j]) for j in range(z.shape[0])])
        base[pr] = names[pc]
    sol = [0] * z.shape[0]
    for i, b in enumerate(base):
        if b in names: sol[names.index(b)] = table[i, -1]
    return sol, z[-1], base

A = [[1, -1], [6, 7], [2, -3]]
b = [-3, 42, 6]
c = [4, -2]
s = ['>=', '<=', '<=']
mode = input("Введіть тип задачі (min/max): ").strip().lower()

A1, b1, s1 = adjust_constraints(A, b, s, mode)
print(f"\nПряма ({mode}):\nF(x) = {expr(c)} → {mode}\nОбмеження:")
for i in range(len(b1)): print(f"  {expr(A1[i])} >= {b1[i]}")

d = dual(A1, b1, c, mode)
print(f"\nДвоїста ({d['type']}):\nZ(y) = {expr(d['c'], 'y')} → {d['type']}\nОбмеження:")
for i in range(len(d['A'])): print(f"  {expr(d['A'][i], 'y')} {d['signs'][i]} {d['b'][i]}")

d['A'], d['b'], d['signs'] = adjust_constraints(d['A'], d['b'], d['signs'], d['type'], flip=True)

tab, z, names, base = simplex_table(d['c'], d['A'], d['b'], d['signs'], d['type'])
sol, opt, _ = solve(tab, z, names, base, d['type'])

print("\nРезультат:")
if sol:
    for i, v in enumerate(sol):
        if v != 0: print(f"y{i+1} = {simplify(v)}")
    print(f"Z(y) = {simplify(opt)}")
else:
    print("Розв’язок не знайдено або необмежений.")
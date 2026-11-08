import sympy as sp
import numpy as np
from sympy import Matrix, simplify, nsimplify

M = sp.Symbol("M")

def form_canonical(c, A, b, signs, goal='max'):
    n = len(c)
    A, b = Matrix(A), Matrix(b).reshape(len(b), 1)
    T = A.row_join(Matrix.eye(len(b)) * 0)
    vars = [f"x{i+1}" for i in range(n)]
    basis, slack, artificial, info = [], 0, 0, []

    for i, s in enumerate(signs):
        if s == "<=":
            v = f"s{slack+1}"
            T[i, n+slack] = 1
            vars.append(v); basis.append(v); slack += 1
        elif s == ">=":
            r, a = f"r{slack+1}", f"a{artificial+1}"
            T[i, n+slack] = -1
            T = T.row_join(Matrix.zeros(T.rows,1)); T[i, T.cols-1] = 1
            vars += [r, a]; basis.append(a)
            slack += 1; artificial += 1; info.append((i,a))
        elif s == "=":
            a = f"a{artificial+1}"
            T = T.row_join(Matrix.zeros(T.rows,1)); T[i, T.cols-1] = 1
            vars.append(a); basis.append(a); artificial += 1; info.append((i,a))

    T = T.row_join(b)
    z = Matrix([sp.Rational(c[i]) if i < len(c) else 0 for i in range(T.cols-1)] + [0])
    mult = 1 if goal == 'min' else -1
    for i, a in info:
        rhs = T[i, -1]
        for j in range(T.cols-1):
            if vars[j] != a: z[j] += mult * M * (-T[i,j])
        z[vars.index(a)] = 0; z[-1] += mult * M * rhs
    z = Matrix([-z[i] if i<len(z)-1 else z[i] for i in range(len(z))])

    return T, z, vars, basis

def identify_pivot_column(z, goal='max'):
    costs = [float(sp.sympify(v).subs(M,1000)) for v in z[:-1]]
    return costs.index(min(costs)) if goal == 'max' and min(costs)<0 else (max([i for i,v in enumerate(costs) if v>0], key=lambda i: costs[i]) if goal == 'min' and any(v>0 for v in costs) else None)

def check_integrality(x): return np.isclose(float(x), round(float(x)))
def extract_fractional_component(x): return float(x) - np.floor(float(x))

def locate_cut_row(T): 
    fracs = [extract_fractional_component(row[-1]) if not check_integrality(row[-1]) else 0 for row in T]
    return (None,0) if all(f==0 for f in fracs) else (np.argmax(fracs), max(fracs))

def cuttin_constraint(T, z, vars, basis):
    r, beta = locate_cut_row(T)
    g_row = [-extract_fractional_component(T[r,j]) for j in range(T.shape[1]-1)] + [1, -extract_fractional_component(T[r,-1])]
    vars.append(f"cut_{len(vars)}"); basis.append(vars[-1])
    T = np.hstack((T[:,:-1], np.zeros((T.shape[0],1)), T[:,-1:]))
    T = np.vstack((T, np.array(g_row, dtype=object)))
    z = np.append(z[:-1],[0,z[-1]])
    return T, z, vars, basis, r, beta, g_row

def execute_SM(T, z, vars, basis, goal='max'):
    T = np.array(T.tolist(), dtype=object); z = [sp.sympify(v) for v in z]
    while (col:=identify_pivot_column(z, goal)) is not None:
        ratios = [T[i,-1]/T[i,col] if T[i,col]>0 else float('inf') for i in range(T.shape[0])]
        if all(r==float('inf') for r in ratios): return None, None, None, None
        row = np.argmin(ratios)
        pivot = T[row,col]
        T[row] = [simplify(x/pivot) for x in T[row]]
        for i in range(T.shape[0]):
            if i!=row:
                T[i] = [simplify(T[i][j] - T[i][col]*T[row][j]) for j in range(len(z))]
        z = [simplify(z[j] - z[col]*T[row][j]) for j in range(len(z))]
        basis[row] = vars[col]
    sol = [0]*(len(z)-1)
    for i in range(T.shape[0]):
        if vars.index(basis[i]) < len(sol):
            sol[vars.index(basis[i])] = T[i,-1]
    return T, z, basis, sol

def execute_dual(T, z, vars, basis):
    while any(T[i,-1]<0 for i in range(T.shape[0])):
        row = np.argmin([T[i,-1] for i in range(T.shape[0])])
        candidates = [(abs(z[j]/T[row,j]), j) for j in range(T.shape[1]-1) if T[row,j]<0]
        if not candidates: return None, None, None, None
        _, col = min(candidates, key=lambda x:float(x[0].subs(M,1000)))
        pivot = T[row,col]
        T[row] = [simplify(x/pivot) for x in T[row]]
        for i in range(T.shape[0]):
            if i!=row:
                T[i] = [simplify(T[i][j]-T[i][col]*T[row][j]) for j in range(T.shape[1])]
        z = [simplify(z[j]-z[col]*T[row][j]) for j in range(len(z))]
        basis[row] = vars[col]
    sol = [0]*(T.shape[1]-1)
    for i in range(T.shape[0]):
        if vars.index(basis[i]) < len(sol):
            sol[vars.index(basis[i])] = T[i,-1]
    return T, z, basis, sol

# === RUN EXAMPLE ===
c = [-2, 1]
A = [[1, 1], [-3, 2]]
b = [6, 3]
signs = ["<=", ">="]
goal = "min"

T, z, vars, basis = form_canonical(c, A, b, signs, goal)
T, z, basis, sol = execute_SM(T, z, vars, basis, goal)

sol_text = ", ".join(f"{vars[i]}={nsimplify(sol[i])}" for i in range(len(c)))
obj_val = nsimplify(sum(sp.Rational(c[i])*sol[i] for i in range(len(c))))

print("1) Simplex solution:\n   ", sol_text, f"\n   F(x) = {obj_val}")

if all(check_integrality(x) for x in sol):
    print("   [✓] Integer solution.")
else:
    print("   [!] Non-integer solution. Applying Gomory cut.")
    T, z, vars, basis, r, beta, g_row = cuttin_constraint(T, z, vars, basis)
    print(f"2) Max fractional part at {basis[r]}: β = {round(beta,4)}")
    print("   Cut coefficients:", [nsimplify(val) for val in g_row])

    T, z, basis, sol = execute_dual(T, z, vars, basis)
    sol_text = ", ".join(f"{vars[i]}={nsimplify(sol[i])}" for i in range(len(c)))
    obj_val = nsimplify(sum(sp.Rational(c[i])*sol[i] for i in range(len(c))))
    print("\n3) Final integer solution:\n   ", sol_text, f"\n   F(x) = {obj_val}")

import sympy as sp
from sympy import Matrix, simplify

M = sp.Symbol("M")
goal = 'max'

numerator = [5, -4]
denominator = [2, 1]
constraints = [
    [5, 2],
    [2, 5],
    [-2, 1],
    [1, 0],
    [0, 1]
]
rhs = [10, 10, 4, 6, 6]
inequalities = [">=", ">=", "<=", "<=", "<="]

def transform_fractional_to_linear(c_num, c_den, mat_constraints, vec_rhs, rel_signs, optimization_goal):
    c_num = list(map(sp.Rational, c_num))
    c_den = list(map(sp.Rational, c_den))
    mat_constraints = [[sp.Rational(entry) for entry in row] for row in mat_constraints]
    vec_rhs = list(map(sp.Rational, vec_rhs))

    var_count = len(c_num)
    linear_objective = c_num + [0]
    new_constraints = []
    new_rhs = []
    new_signs = []

    for idx in range(len(mat_constraints)):
        if vec_rhs[idx] < 0:
            mat_constraints[idx] = [-val for val in mat_constraints[idx]]
            vec_rhs[idx] = -vec_rhs[idx]
            rel_signs[idx] = ">=" if rel_signs[idx] == "<=" else "<=" if rel_signs[idx] == ">=" else rel_signs[idx]

        constraint_row = mat_constraints[idx] + [-vec_rhs[idx]]

        if rel_signs[idx] == ">=":
            constraint_row = [-val for val in constraint_row]
            new_signs.append("<=")
        else:
            new_signs.append(rel_signs[idx])

        new_constraints.append(constraint_row)
        new_rhs.append(0)

    new_constraints.append(c_den + [0])
    new_signs.append("=")
    new_rhs.append(1)

    var_names = [f"y{i+1}" for i in range(var_count)] + ["y0"]
    return new_constraints, linear_objective, var_names, new_signs, new_rhs

def display_fractional_objective(c_num, c_den):
    def format_expr(coeffs):
        terms = []
        for i, coeff in enumerate(coeffs):
            var = f"x{i+1}"
            if coeff > 0 and terms:
                terms.append(f"+ {coeff}{var}")
            elif coeff < 0:
                terms.append(f"- {abs(coeff)}{var}")
            else:
                terms.append(f"{coeff}{var}")
        return " ".join(terms)

    num_expr = format_expr(c_num)
    den_expr = format_expr(c_den)

    print("Дробово-лінійна цільова функція:")
    print(f"F(x) = ({num_expr}) / ({den_expr})")

def display_linearized_problem(obj, A, b, variables, signs, goal):
    sym_vars = [sp.Symbol(v) for v in variables]
    expr = sp.Add(*[a * v for a, v in zip(obj, sym_vars)])
    print(f"\nЛінеаризована задача:\n{goal}(G) = {expr}")
    for i in range(len(A)):
        lhs = sp.Add(*[A[i][j]*sym_vars[j] for j in range(len(variables))])
        print(f"   {lhs} {signs[i]} {b[i]}")

def create_tableau(obj, A, b, signs, goal):
    A = Matrix(A)
    b = Matrix(b).reshape(len(b), 1)
    T = A
    names = [f"y{i+1}" for i in range(A.shape[1] - 1)] + ["y0"]
    base, artificial = [], []
    s_id = a_id = 1

    for i, s in enumerate(signs):
        if s == "<=":
            T = T.row_join(Matrix.eye(len(signs)).col(i))
            var = f"s{s_id}"; names.append(var); base.append(var)
            s_id += 1
        elif s == ">=":
            T = T.row_join(-Matrix.eye(len(signs)).col(i))
            T = T.row_join(Matrix.eye(len(signs)).col(i))
            names += [f"s{s_id}", f"w{a_id}"]
            base.append(f"w{a_id}")
            artificial.append((i, f"w{a_id}"))
            s_id += 1; a_id += 1
        elif s == "=":
            T = T.row_join(Matrix.eye(len(signs)).col(i))
            var = f"w{a_id}"
            names.append(var); base.append(var); artificial.append((i, var))
            a_id += 1

    T = T.row_join(b)
    c = Matrix([sp.Rational(obj[i]) if i < len(obj) else 0 for i in range(T.cols)])

    mult = 1 if goal == 'min' else -1
    for idx, av in artificial:
        row = T[idx, :]
        for j in range(len(c) - 1):
            if names[j] != av:
                c[j] += mult * M * (-row[j])
        c[names.index(av)] = 0
        c[-1] += mult * M * row[-1]

    for i in range(len(c) - 1): c[i] *= -1
    return T, c, names, base

def solve_with_simplex(T, c, names, base, goal):
    while True:
        row = [float(v.subs(M, 1000)) for v in c[:-1]]
        col = (
            row.index(min(row)) if min(row) < 0 else None
        ) if goal == 'max' else (
            max((i for i, v in enumerate(row) if v > 0), key=lambda i: row[i], default=None)
        )
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

    solution = [0] * c.shape[0]
    for i in range(T.shape[0]):
        solution[names.index(base[i])] = T[i, -1]
    return solution, c[-1], base

display_fractional_objective(numerator, denominator)
A_lin, c_lin, vars_lin, signs_lin, b_lin = transform_fractional_to_linear(numerator, denominator, constraints, rhs, inequalities, goal)
display_linearized_problem(c_lin, A_lin, b_lin, vars_lin, signs_lin, goal)
T, c_row, names, base = create_tableau(c_lin, A_lin, b_lin, signs_lin, goal)
solution, final_val, _ = solve_with_simplex(T.copy(), c_row.copy(), names, base, goal)

if solution:
    y0_val = solution[names.index("y0")]

    print("\nЗначення змінних y:")
    y_values = "; ".join([f"{var} = {simplify(solution[names.index(var)])}" for var in ["y1", "y2", "y0"]])
    print(f"   {y_values}")

    if y0_val == 0:
        print("\ny0 тотожний нулю, зворотня заміна неможлива.")
    else:
        print("\nЗворотня заміна ( xi = yi / y0):")
        xs = [simplify(solution[names.index(f"y{i+1}")] / y0_val) for i in range(len(numerator))]
        x_values = "; ".join([f"x{i+1} = {xi}" for i, xi in enumerate(xs)])
        print(f"   {x_values}")

        f_val = simplify(sum(numerator[i] * xs[i] for i in range(len(xs))) / sum(denominator[i] * xs[i] for i in range(len(xs))))
        print(f"\nЗначення дробово-лінійної функції: F(x) = {f_val}")
else:
    print("\nРозв’язок не знайдено.")
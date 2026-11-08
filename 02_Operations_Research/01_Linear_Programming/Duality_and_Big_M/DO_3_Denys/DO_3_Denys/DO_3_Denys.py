import numpy as np

M = 10**6

def format_expression(coeffs, prefix='x'):
    return " ".join(f"{'+' if coef >= 0 and i > 0 else ''}{coef}{prefix}{i+1}" for i, coef in enumerate(coeffs))

def build_dual(A, b, c, sense):
    A = np.array(A)
    return {
        'dual_sense': 'max' if sense == 'min' else 'min',
        'A_dual': A.T.tolist(),
        'b_dual': c,
        'c_dual': b,
        'constraint_signs': ['<=' if sense == 'min' else '>='] * A.shape[1]
    }

def print_problem(A, b, c, sense):
    print(f"{'='*25} {sense.upper()} {'='*25}")
    print(f"F = {format_expression(c)} {'→ min' if sense == 'min' else '→ max'}")
    print("Система обмежень:")
    for row, val in zip(A, b):
        print(f"  {format_expression(row)} = {val}")

def print_dual(dual):
    print(f"\nZ = {format_expression(dual['c_dual'], 'y')} {'→ min' if dual['dual_sense'] == 'min' else '→ max'}")
    print("Система обмежень:")
    for row, val, sign in zip(dual['A_dual'], dual['b_dual'], dual['constraint_signs']):
        print(f"  {format_expression(row, 'y')} {sign} {val}")

def initialize_table(goal):
    return np.array([
        [1, 3, -1, 0, 1, 4] if goal == 'max' else [-1, -3, 1, 0, 0, 4],
        [-1, 2, 0, 1, 0, 2] if goal == 'max' else [1, -2, 0, -1, 1, 2]
    ], dtype=float)

def build_objective(goal):
    return [-M + 6, -3*M - 3, M, 0, 0, -4*M] if goal == 'max' else [M - 6, -2*M + 3, 0, -M, 0, 2*M]

def choose_pivot(c, goal):
    return np.argmin(c[:-1]) if goal == 'max' and np.min(c[:-1]) < 0 else (
        max((j for j in range(len(c)-1) if c[j] > 0), key=lambda j: c[j], default=None) if goal == 'min' else None
    )

def simplex(table, c, goal):
    num_vars, num_constraints = len(c) - 1, len(table)
    basis = [col for col in range(num_vars) if list(table[:, col]).count(1) == 1 and list(table[:, col]).count(0) == num_constraints - 1][:num_constraints]

    while (pivot := choose_pivot(c, goal)) is not None:
        ratios = [table[i, -1] / table[i, pivot] if table[i, pivot] > 0 else np.inf for i in range(num_constraints)]
        if all(r == np.inf for r in ratios): return None, None, None
        row = np.argmin(ratios)
        basis[row] = pivot
        table[row] /= table[row, pivot]
        for i in range(num_constraints):
            if i != row:
                table[i] -= table[i, pivot] * table[row]
        c -= c[pivot] * table[row]

    sol = np.zeros(num_vars)
    for i in range(num_constraints):
        sol[basis[i]] = table[i, -1]
    return sol, c[-1], basis

def print_solution(solution, value):
    print("\nОптимальне рішення:")
    print("(" + ", ".join(f"y{i+1} = {solution[i]:.2f}" for i in range(len(solution))) + ")")
    print(f"Значення цільової функції при цьому плані: {value:.2f}")

examples = [
    {'A': [[-1, -1], [-3, 2]], 'b': [-6, 3], 'c': [-4, 2], 'sense': 'min'},
    {'A': [[1, 1], [3, -2]], 'b': [6, -3], 'c': [-4, 2], 'sense': 'max'}
]

for ex in examples:
    print_problem(ex['A'], ex['b'], ex['c'], ex['sense'])
    dual = build_dual(ex['A'], ex['b'], ex['c'], ex['sense'])
    print_dual(dual)
    goal = dual['dual_sense']
    table = initialize_table(goal)
    obj_func = build_objective(goal)
    sol, val, _ = simplex(table.copy(), obj_func, goal)
    print_solution(sol, val)
import numpy as np

def build_dual(A, b, c, sense='min'):
    A = np.array(A)
    b = np.array(b)
    c = np.array(c)
    A_dual = A.T
    b_dual = c
    c_dual = b
    dual_sense = 'max' if sense == 'min' else 'min'
    constraint_signs = ['<='] * A_dual.shape[0] if sense == 'min' else ['>='] * A_dual.shape[0]
    return {
        'dual_sense': dual_sense,
        'A_dual': A_dual.tolist(),
        'b_dual': b_dual.tolist(),
        'c_dual': c_dual.tolist(),
        'constraint_signs': constraint_signs
    }

def format_expression(coeffs, vars_prefix='x'):
    return " ".join([
        f"{'+' if coef >= 0 and i > 0 else ''}{coef}{vars_prefix}{i+1}"
        for i, coef in enumerate(coeffs)
    ])

def format_lp_problem(A, b, c, sense, name):
    sense_symbol = "→ min" if sense == 'min' else "→ max"
    output = [f"🔹 Пряма задача: {name}",
              f"F = {format_expression(c)} {sense_symbol}"]
    output.append("За умов:")
    for i, row in enumerate(A):
        output.append(f"  {format_expression(row)} {b[i]}")
    return "\n".join(output)

def format_dual_output(dual):
    sense_symbol = "→ min" if dual['dual_sense'] == 'min' else "→ max"
    output = [f"🔸 Двоїста задача",
              f"Z = {format_expression(dual['c_dual'], vars_prefix='y')} {sense_symbol}"]
    output.append("За умов:")
    for i in range(len(dual['A_dual'])):
        row_expr = format_expression(dual['A_dual'][i], vars_prefix='y')
        output.append(f"  {row_expr} {dual['constraint_signs'][i]} {dual['b_dual'][i]}")
    return "\n".join(output)

# ================= ВХІДНІ ДАНІ =================

A_min = [[-2, 1], [-1, -1], [-3, 2]]
b_min = [-8, -5, 3]
c_min = [-4, 2]
dual_min = build_dual(A_min, b_min, c_min, sense='min')

A_max = [[2, -1], [1, 1], [3, -2]]
b_max = [8, 5, -3]
c_max = [-4, 2]
dual_max = build_dual(A_max, b_max, c_max, sense='max')

# ================= ВИВІД =================

print(format_lp_problem(A_min, b_min, c_min, "min", "Приклад 1"))
print()
print(format_dual_output(dual_min))
print("\n" + "="*50 + "\n")

print(format_lp_problem(A_max, b_max, c_max, "max", "Приклад 2"))
print()
print(format_dual_output(dual_max))
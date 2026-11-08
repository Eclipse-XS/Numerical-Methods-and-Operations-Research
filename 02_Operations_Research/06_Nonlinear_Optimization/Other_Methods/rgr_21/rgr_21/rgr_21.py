import numpy as np
from sympy import symbols, lambdify, sqrt

# Символьні змінні
x1, x2, lam = symbols('x1 x2 λ')

# Функція
f_expr = 8*x1**2 - 6*x1*x2 + 3*x2**2 - 7*x1 + 10*x2

# Градієнт
df_dx1 = f_expr.diff(x1)
df_dx2 = f_expr.diff(x2)

# Lambdify
f_lambdified = lambdify((x1, x2), f_expr, 'numpy')
grad_f = lambdify((x1, x2), [df_dx1, df_dx2], 'numpy')

# Функція вздовж напрямку
def line_function(x_k, d_k):
    x_lambda = x_k + lam * d_k
    f_lambda_expr = f_expr.subs({x1: x_lambda[0], x2: x_lambda[1]})
    return f_lambda_expr.expand()

# Знаходження мінімуму параболи
def find_parabola_min(f_lambda_expr):
    a = f_lambda_expr.coeff(lam, 2)
    b = f_lambda_expr.coeff(lam, 1)
    return float(b) / (2 * float(a))

# Початкові умови
x_k = np.array([1.5, 2.5])
eps = 0.01

for k in range(10):  # максимум 100 ітерацій
    grad_k = np.array(grad_f(x_k[0], x_k[1]))
    norm_grad = np.linalg.norm(grad_k)

    print(f"ІТЕРАЦІЯ №{k+1}")
    print(f"x⁽{k}⁾ = ({x_k[0]:.4f}, {x_k[1]:.4f})")
    print(f"∇f(x⁽{k}⁾) = ({grad_k[0]:.4f}, {grad_k[1]:.4f})")
    print(f"||∇f|| = {norm_grad:.4f}")

    if norm_grad < eps:
        break

    d_k = -grad_k
    print(f"d = ({d_k[0]:.4f}, {d_k[1]:.4f})")

    f_lambda_expr = line_function(x_k, d_k)
    lam_k = find_parabola_min(f_lambda_expr)
    print(f"λ = {lam_k:.4f}")

    x_k = x_k + lam_k * d_k
    print(f"x⁽{k+1}⁾ = ({x_k[0]:.4f}, {x_k[1]:.4f})")
    print("-" * 50)

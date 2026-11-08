import numpy as np
import sympy as sp

# Use sympy symbols for full symbolic trace
x1, x2, alpha = sp.symbols('x1 x2 alpha')

# Define the function symbolically
f_expr = 2*x1**2 - 7*x1*x2 + x2**2 - 2*x1 + 3*x2

# Compute gradient
grad_f_expr = [sp.diff(f_expr, x1), sp.diff(f_expr, x2)]

# Function to compute numeric gradient at a given point
def compute_gradient(x_val):
    subs = {x1: x_val[0], x2: x_val[1]}
    grad_numeric = np.array([float(expr.evalf(subs=subs)) for expr in grad_f_expr])
    print(f"Gradient at point {x_val}:")
    print(f" ∂f/∂x1 = {grad_f_expr[0]} → {grad_numeric[0]}")
    print(f" ∂f/∂x2 = {grad_f_expr[1]} → {grad_numeric[1]}")
    return grad_numeric

# Function to perform line minimization along direction d
def line_minimization(f_expr, x_val, d):
    # substitute x1 = x1_0 + alpha*d1, x2 = x2_0 + alpha*d2 into f
    f_alpha = f_expr.subs({x1: x_val[0] + alpha * d[0], x2: x_val[1] + alpha * d[1]})
    f_alpha_simplified = sp.simplify(sp.expand(f_alpha))
    print(f"\nFunction along direction d: f(α) = {f_alpha_simplified}")
    # take derivative and solve
    f_prime = sp.diff(f_alpha_simplified, alpha)
    alpha_opt = float(sp.solve(f_prime, alpha)[0])
    print(f" → df/dα = {f_prime}")
    print(f" → Optimal α = {alpha_opt}")
    return alpha_opt

# Perform one iteration of DFP
def dfp_iteration(x0, A0):
    print(f"\n--- Starting DFP iteration ---")
    print(f"x⁽⁰⁾ = {x0}")
    g0 = compute_gradient(x0)
    d0 = -A0 @ g0
    print(f"\nDirection d⁽⁰⁾ = -A⁽⁰⁾ * grad = {d0}")

    alpha0 = line_minimization(f_expr, x0, d0)
    x1 = x0 + alpha0 * d0
    print(f"\nNew point x⁽¹⁾ = x⁽⁰⁾ + α * d = {x1}")

    g1 = compute_gradient(x1)
    v = x1 - x0
    u = g1 - g0

    print(f"\nΔx = {v}")
    print(f"Δg = {u}")

    # DFP matrix update
    A1 = A0 + np.outer(v, v) / (v @ u) - (A0 @ np.outer(u, u) @ A0) / (u @ A0 @ u)
    print(f"\nUpdated matrix A⁽¹⁾:\n{A1}")

    return x1, A1, g1

# Initial conditions
x0 = np.array([0.1, 0.1])
A0 = np.eye(2)

# Run one DFP iteration
x1, A1, g1 = dfp_iteration(x0, A0)


import numpy as np
import matplotlib.pyplot as plt

h = 0.25
x_start, x_end = 0, 1
n = int((x_end - x_start) / h) + 1

x = np.linspace(x_start, x_end, n)

def p(x):
    return x

def q(x):
    return -1

def f(x):
    return x**2 - x + 1

A = np.zeros((n, n))
b = np.zeros(n)

A[0, 0] = (h - 1) / h
A[0, 1] = 1 / h
b[0] = 1.5

for k in range(1, n-1):
    A[k, k-1] = 1 - (h / 2) * p(x[k])
    A[k, k] = h**2 * q(x[k]) - 2
    A[k, k+1] = 1 + (h / 2) * p(x[k])
    b[k] = h**2 * f(x[k])

A[-1, -2] = -1 / h
A[-1, -1] = 1 / h
b[-1] = 1

for i in range(n):
    if A[i, i] == 0:
        A[i, i] = 1e-10

y = np.linalg.solve(A, b)

for i in range(n):
    print(f"x = {x[i]:.2f}, y = {y[i]:.5f}")

plt.plot(x, y, marker='o', linestyle='-', color='b')
plt.xlabel("x")
plt.ylabel("y")
plt.title("Розв'язок диференціального рівняння методом сіток")
plt.grid(True)
plt.show()

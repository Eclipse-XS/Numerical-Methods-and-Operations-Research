import numpy as np

# 🔹 Цільова функція
def f(x):
    x1, x2 = x
    return 4 * x1**2 - 4 * x1 * x2 + 3 * x2**2 + x1

# 🔹 Градієнт функції
def grad_f(x):
    x1, x2 = x
    df_dx1 = 8 * x1 - 4 * x2 + 1
    df_dx2 = -4 * x1 + 6 * x2
    return np.array([df_dx1, df_dx2])

# 🔹 Метод Флетчера–Рівса
def fletcher_reeves(x0, epsilon=1e-3, max_iter=100):
    x = x0
    grad = grad_f(x)
    d = -grad
    history = [(x.copy(), f(x), grad.copy())]

    for k in range(max_iter):
        # Одновимірна мінімізація вздовж напрямку d
        def phi(t):
            return f(x + t * d)

        # Пошук t_best на відрізку [0, 1]
        t_vals = np.linspace(0, 1, 1000)
        t_best = t_vals[np.argmin([phi(t) for t in t_vals])]

        # Оновлення точки
        x_new = x + t_best * d
        grad_new = grad_f(x_new)

        # Зупинка
        if np.linalg.norm(grad_new) < epsilon:
            history.append((x_new.copy(), f(x_new), grad_new.copy()))
            break

        beta = np.dot(grad_new, grad_new) / np.dot(grad, grad)
        d = -grad_new + beta * d

        x = x_new
        grad = grad_new
        history.append((x.copy(), f(x), grad.copy()))

    return history

# 🔹 Запуск алгоритму
if __name__ == "__main__":
    x0 = np.array([0.0, 0.0])
    result = fletcher_reeves(x0)

    print("🔁 Ітерації методу Флетчера–Рівса:\n")
    for i, (x, fx, g) in enumerate(result):
        norm_g = np.linalg.norm(g)
        print(f"Ітерація {i}:")
        print(f"  x = ({x[0]:.6f}, {x[1]:.6f})")
        print(f"  f(x) = {fx:.6f}")
        print(f"  ∇f(x) = ({g[0]:.6f}, {g[1]:.6f}), ‖∇f‖ = {norm_g:.6f}\n")

    print("✅ Глобальний мінімум наближено знайдено.")


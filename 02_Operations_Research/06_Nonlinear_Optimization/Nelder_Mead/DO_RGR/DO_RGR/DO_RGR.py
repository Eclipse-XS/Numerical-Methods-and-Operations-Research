import numpy as np

# 🔹 Цільова функція
def f(x, y):
    return -x**2 - x*y - 3*y**2 + 12*x + 11*y - 6

# 🔹 Формули δ1 і δ2 для регулярного симплексу
def compute_deltas(alpha):
    delta1 = ((np.sqrt(3) + 1) / (2 * np.sqrt(2))) * alpha
    delta2 = ((np.sqrt(3) - 1) / (2 * np.sqrt(2))) * alpha
    return delta1, delta2

# 🔹 Початковий симплекс (у 2D — три точки)
def initial_simplex(x0, alpha):
    delta1, delta2 = compute_deltas(alpha)
    x1 = x0 + np.array([delta2, delta1])
    x2 = x0 + np.array([delta1, delta2])
    return [np.array(x0, dtype=float), x1, x2]

# 🔹 Центр тяжіння
def centroid(points):
    return np.mean(points, axis=0)

# 🔹 Метод симплекс-пошуку
def simplex_search(f, x0, alpha=2, max_iter=10):
    simplex = initial_simplex(np.array(x0, dtype=float), alpha)

    for iteration in range(max_iter):
        print(f"\n🔁 Ітерація {iteration + 1}")

        # Обчислення значень функції
        values = [f(p[0], p[1]) for p in simplex]

        # Вивід симплексу
        for i, point in enumerate(simplex):
            print(f"  Точка {i}: {point}, f = {values[i]:.5f}")

        # Сортування (по спаданню)
        sorted_points = sorted(zip(simplex, values), key=lambda x: -x[1])
        best = sorted_points[0][0]
        middle = sorted_points[1][0]
        worst = sorted_points[2][0]

        print(f"  🔹 Найкраща: {best}, f = {f(best[0], best[1]):.5f}")
        print(f"  🔻 Найгірша: {worst}, f = {f(worst[0], worst[1]):.5f}")

        # Центр тяжіння
        x_c = centroid([best, middle])
        print(f"  🧲 Центр тяжіння: {x_c}")

        # Відображення
        x_r = 2 * x_c - worst
        f_r = f(x_r[0], x_r[1])
        print(f"  📍 Відображення: {x_r}, f = {f_r:.5f}")

        # Оновлення симплексу
        if f_r > f(worst[0], worst[1]):
            print("  ✅ Відображення прийнято — замінюємо найгіршу точку")
            for i, point in enumerate(simplex):
                if np.allclose(point, worst):
                    simplex[i] = x_r
                    break
        else:
            print("  ⚠️  Відображення гірше — стискаємо симплекс навколо найкращої точки")
            simplex = [best + 0.5 * (p - best) for p in simplex]

    # Після завершення
    best = max(simplex, key=lambda p: f(p[0], p[1]))
    return best, f(best[0], best[1])

# 🔹 Запуск
if __name__ == "__main__":
    best_point, best_value = simplex_search(f, x0=[0, 0], alpha=2, max_iter=10)
    print("\n🔹 Оптимум знайдено:")
    print(f"   x = {best_point[0]:.4f}, y = {best_point[1]:.4f}")
    print(f"   f(x, y) = {best_value:.5f}")

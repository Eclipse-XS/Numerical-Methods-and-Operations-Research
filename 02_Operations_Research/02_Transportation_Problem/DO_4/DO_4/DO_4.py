import numpy as np
from collections import Counter

def balance_problem(supply, demand, costs):
    total_supply = sum(supply)
    total_demand = sum(demand)

    if total_supply == total_demand:
        print("Задача є закритою (збалансованою).")
        return supply, demand, costs, False, None

    if total_supply < total_demand:
        diff = total_demand - total_supply
        print(f"Задача є відкритою. Додається фіктивний постачальник з обсягом {diff}.")
        supply = np.append(supply, diff)
        fake_row = np.zeros((1, costs.shape[1]))
        costs = np.vstack([costs, fake_row])
        return supply, demand, costs, True, 'supplier'
    else:
        diff = total_supply - total_demand
        print(f"Задача є відкритою. Додається фіктивний споживач з потребою {diff}.")
        demand = np.append(demand, diff)
        fake_col = np.zeros((costs.shape[0], 1))
        costs = np.hstack([costs, fake_col])
        return supply, demand, costs, True, 'consumer'

def northwest_corner_method(a, b):
    a = a.copy()
    b = b.copy()
    m, n = len(a), len(b)
    x = np.zeros((m, n), dtype=int)
    i = j = 0
    while i < m and j < n:
        value = min(a[i], b[j])
        x[i, j] = value
        a[i] -= value
        b[j] -= value
        if a[i] == 0 and b[j] == 0:
            if j + 1 < n:
                j += 1
            elif i + 1 < m:
                i += 1
            else:
                break
        elif a[i] == 0:
            i += 1
        elif b[j] == 0:
            j += 1
    return x

def least_cost_method(C, a, b):
    a = a.copy()
    b = b.copy()
    m, n = len(a), len(b)
    x = np.zeros((m, n), dtype=int)
    cost_indices = [(i, j) for i in range(m) for j in range(n)]
    cost_indices.sort(key=lambda x: C[x[0]][x[1]])
    for i, j in cost_indices:
        if a[i] == 0 or b[j] == 0:
            continue
        value = min(a[i], b[j])
        x[i, j] = value
        a[i] -= value
        b[j] -= value
    return x

def transport(supply, demand, costs, initial_plan):
    C = np.copy(costs)
    X = initial_plan.astype(float)
    n, m = C.shape

    _x, _y = np.where(X > 0)
    nonzero = list(zip(_x, _y))
    if len(nonzero) < n + m - 1:
        for i in range(n):
            for j in range(m):
                if X[i, j] == 0:
                    X[i, j] = 0.00001
                    nonzero.append((i, j))
                    if len(nonzero) >= n + m - 1:
                        break
            if len(nonzero) >= n + m - 1:
                break

    while True:
        u = np.array([np.nan] * n)
        v = np.array([np.nan] * m)
        S = np.zeros((n, m))

        _x, _y = np.where(X > 0)
        nonzero = list(zip(_x, _y))
        u[nonzero[0][0]] = 0

        changed = True
        while changed:
            changed = False
            for i, j in nonzero:
                if np.isnan(u[i]) and not np.isnan(v[j]):
                    u[i] = C[i, j] - v[j]
                    changed = True
                elif not np.isnan(u[i]) and np.isnan(v[j]):
                    v[j] = C[i, j] - u[i]
                    changed = True

        for i in range(n):
            for j in range(m):
                S[i, j] = C[i, j] - u[i] - v[j]

        if np.all(S >= 0):
            break

        i, j = np.unravel_index(np.argmin(S), S.shape)
        start = (i, j)
        T = np.copy(X)
        T[start] = 1

        max_iters = 1000
        iter_count = 0
        while True:
            iter_count += 1
            if iter_count > max_iters:
                return X, np.sum(X * C)

            _xs, _ys = np.nonzero(T)
            xcount, ycount = Counter(_xs), Counter(_ys)
            for x, count in xcount.items():
                if count <= 1:
                    T[x, :] = 0
            for y, count in ycount.items():
                if count <= 1:
                    T[:, y] = 0
            if all(count > 1 for count in xcount.values()) and all(count > 1 for count in ycount.values()):
                break

        dist = lambda a, b: abs(a[0] - b[0]) + abs(a[1] - b[1])
        fringe = set(tuple(p) for p in np.argwhere(T > 0))
        size = len(fringe)
        path = [start]
        while len(path) < size:
            last = path[-1]
            if last in fringe:
                fringe.remove(last)
            next_candidates = [f for f in fringe if f[0] == last[0] or f[1] == last[1]]
            if not next_candidates:
                return X, np.sum(X * C)
            next_step = min(next_candidates, key=lambda x: dist(last, x))
            path.append(next_step)

        neg = path[1::2]
        pos = path[::2]

        if not neg:
            return X, np.sum(X * C)

        q = min([X[i, j] for i, j in neg])
        for i, j in neg:
            X[i, j] -= q
        for i, j in pos:
            X[i, j] += q

    return X, np.sum(X * C)

if __name__ == '__main__':
    supply = np.array([230, 90, 40])
    demand = np.array([110, 30, 50, 80, 90])
    costs = np.array([[4, 5, 6, 8, 10], [10, 3, 2, 3, 15], [4, 10, 5, 1, 16]])

    supply, demand, costs, is_open, added_type = balance_problem(supply, demand, costs)

    print("\nОберіть метод для побудови початкового плану:")
    print("1 - Північно-західний кут")
    print("2 - Метод найменшої вартості")
    method = input("Ваш вибір (1 або 2): ")

    if method == '1':
        initial_plan = northwest_corner_method(supply, demand)
        method_name = "Північно-західний кут"
    elif method == '2':
        initial_plan = least_cost_method(costs, supply, demand)
        method_name = "Метод найменшої вартості"
    else:
        print("Некоректний вибір методу.")
        exit()

    print(f"\nПочатковий опорний план ({method_name}):")
    print(initial_plan)
    initial_cost = np.sum(initial_plan * costs)
    print(f"Початкова вартість перевезення: {initial_cost}")

    plan, total_cost = transport(supply, demand, costs, initial_plan)

    if plan is not None:
        print("\nОптимізований план перевезень:")
        plan_int = np.round(plan).astype(int)

        if is_open and added_type == 'supplier':
            plan_int = plan_int[:-1, :]
        elif is_open and added_type == 'consumer':
            plan_int = plan_int[:, :-1]

        print(plan_int)
        print(f"Мінімальна вартість перевезення: {int(round(total_cost))}")

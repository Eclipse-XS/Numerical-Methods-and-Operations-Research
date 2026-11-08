import numpy as np
from collections import Counter
from tabulate import tabulate
import os

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')

def show_transport_table(costs, supply, demand, plan=None, title="Transport Table"):
    headers = [f"B{j+1}" for j in range(len(demand))] + ["Supply"]
    table = []

    for i in range(len(supply)):
        row = [f"A{i+1}"]
        for j in range(len(demand)):
            if plan is not None and plan[i, j] > 0:
                row.append(f"{int(costs[i, j])} [{int(plan[i, j])}]")
            else:
                row.append("-" if costs[i, j] == 0 else str(int(costs[i, j])))
        row.append(int(supply[i]))
        table.append(row)

    footer = ["Demand"] + list(map(int, demand)) + [""]
    print(f"\n{title}")
    print(tabulate(table + [footer], headers=[""] + headers, tablefmt="fancy_grid"))

def validate_balance(supply, demand, costs):
    total_supply = sum(supply)
    total_demand = sum(demand)

    if total_supply == total_demand:
        print("- Задача закрита.")
        return supply, demand, costs, False, None

    diff = abs(total_supply - total_demand)
    print("- Задача відкрита.")

    if total_supply < total_demand:
        print(f"- Додаємо фіктивного постачальника")
        supply = np.append(supply, diff)
        costs = np.vstack([costs, np.zeros((1, costs.shape[1]))])
        return supply, demand, costs, True, 'supplier'
    else:
        print(f"- Додаємо фіктивного споживача")
        demand = np.append(demand, diff)
        costs = np.hstack([costs, np.zeros((costs.shape[0], 1))])
        return supply, demand, costs, True, 'consumer'

def northwest(stock, need):
    stock = stock.copy()
    need = need.copy()
    suppliers = len(stock)
    consumers = len(need)

    plan = np.full((suppliers, consumers), np.nan)
    i = j = 0

    while i < suppliers and j < consumers:
        shipment = min(stock[i], need[j])
        plan[i, j] = shipment
        stock[i] -= shipment
        need[j] -= shipment

        if stock[i] == 0 and need[j] == 0:
            if j + 1 < consumers:
                j += 1
            elif i + 1 < suppliers:
                i += 1
            else:
                break
        elif stock[i] == 0:
            i += 1
        elif need[j] == 0:
            j += 1

    return np.nan_to_num(plan, nan=0).astype(int)

def cheapest_value(C, a, b):
    a, b = a.copy(), b.copy()
    m, n = len(a), len(b)
    x = np.zeros((m, n), dtype=int)

    pos = [(i, j) for i in range(m) for j in range(n) if C[i, j] > 0]
    zero = [(i, j) for i in range(m) for j in range(n) if C[i, j] == 0]

    pos.sort(key=lambda p: (C[p], p[0], p[1]))
    zero.sort(key=lambda p: (p[0], p[1]))

    for group in [pos, zero]:
        for i, j in group:
            if a[i] == 0 or b[j] == 0:
                continue
            val = min(a[i], b[j])
            x[i, j] = val
            a[i] -= val
            b[j] -= val

    return x

def potentials(supply, demand, costs, plan):
    C = np.copy(costs)
    X = plan.astype(float)
    n, m = C.shape

    while True:
        u = np.array([np.nan] * n)
        v = np.array([np.nan] * m)
        S = np.zeros((n, m))

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
        q = min([X[i, j] for i, j in neg])
        for i, j in neg:
            X[i, j] -= q
        for i, j in pos:
            X[i, j] += q

    return X, np.sum(X * C)

def calculate_cost(c, x):
    return int(sum(c[i, j] * x[i, j] for i in range(c.shape[0]) for j in range(c.shape[1])))

def lab_4():
    print("Оберіть метод побудови початкового плану транспортування:")
    print("1 — Правило північно-західного кута")
    print("2 — Метод найменших витрат")
    selected = input("Ваш вибір (1 або 2): ").strip()

    if selected == "1":
        method_code = "NWC"
    elif selected == "2":
        method_code = "LC"
    else:
        print("Невірний вибір. Завершення програми.")
        return

    clear_console()

    stock = np.array([30, 40, 70, 60])
    need = np.array([30, 80, 20, 70])

    tariffs = np.array([
        [1, 9, 7, 2],
        [3, 1, 5, 5],
        [6, 8, 3, 4],
        [2, 3, 1, 3]
    ])

    print("\n--- КРОК 1: Перевірка збалансованості -------------------------")
    stock, need, tariffs, *_ = validate_balance(stock, need, tariffs)
    if method_code == "NWC":
        initial_plan = northwest(stock, need)
        method_title = "Правило північно-західного кута"
    elif method_code == "LC":
        initial_plan = cheapest_value(tariffs, stock, need)
        method_title = "Метод найменших витрат"

    print(f"\n--- КРОК 2: Початковий план ({method_title}) ----------------------------")
    show_transport_table(tariffs, stock, need, plan=initial_plan)
    initial_cost = calculate_cost(tariffs, initial_plan)
    print(f"\nF(x) (початковий план): {initial_cost}")

    print("\n--- КРОК 3: Оптимізація методом потенціалів -------------------")
    optimized_plan, optimized_cost = potentials(stock, need, tariffs, initial_plan)
    show_transport_table(
        tariffs[:len(stock), :len(need)],
        stock,
        need,
        plan=np.round(optimized_plan).astype(int)
    )
    print(f"\nF(x) (оптимізований план): {int(optimized_cost)}")

    print("\n--- ПОРІВНЯННЯ РЕЗУЛЬТАТІВ -------------------------------------")
    print(f"F(x) початкового плану       = {initial_cost}")
    print(f"F(x) після оптимізації       = {int(optimized_cost)}")

if __name__ == '__main__':
    lab_4()
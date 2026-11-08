import numpy as np
from collections import Counter
from tabulate import tabulate

def print_transport_table(costs, supply, demand, title="Транспортна таблиця", plan=None):
    headers = [f"B{j+1}" for j in range(len(demand))] + ["aᵢ"]
    table = []

    for i in range(len(supply)):
        row = []
        for j in range(len(demand)):
            cost = int(costs[i][j])
            if plan is not None and plan[i][j] > 0:
                row.append(f"{cost}[{int(plan[i][j])}]")
            else:
                row.append("-" if cost == 0 else str(cost))
        row.append(int(supply[i]))
        table.append([f"A{i+1}"] + row)

    b_row = ["bⱼ"] + [int(val) for val in demand] + [""]
    print(f"\n--- {title} ---")
    print(tabulate(table + [b_row], headers=[""] + headers, tablefmt="fancy_grid"))

def check_problem_type(supply, demand, costs):
    total_supply = sum(supply)
    total_demand = sum(demand)

    print("\n======= ТРАНСПОРТНА ЗАДАЧА =======")
    print("=== Балансування задачі ===")

    if total_supply == total_demand:
        print("Задача є закритою (збалансованою).")
        return supply, demand, costs, False, None

    if total_supply < total_demand:
        diff = total_demand - total_supply
        print(f"До ТЗ-таблиці додаємо фіктивного постачальника із запасом ai={diff}.")
        supply = np.append(supply, diff)
        fake_row = np.zeros((1, costs.shape[1]))
        costs = np.vstack([costs, fake_row])
        return supply, demand, costs, True, 'supplier'
    else:
        diff = total_supply - total_demand
        print(f"До ТЗ-таблиці додаємо фіктивного споживача із потребою bj={diff}.")
        demand = np.append(demand, diff)
        fake_col = np.zeros((costs.shape[0], 1))
        costs = np.hstack([costs, fake_col])
        return supply, demand, costs, True, 'consumer'

def northwest_corner_method(a, b):
    a = a.copy()
    b = b.copy()
    m, n = len(a), len(b)
    x = np.full((m, n), np.nan)
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

    x = np.nan_to_num(x, nan=0).astype(int)
    return x

def least_cost_method(C, a, b):
    a = a.copy()
    b = b.copy()
    m, n = len(a), len(b)
    x = np.zeros((m, n), dtype=int)

    positive_costs = [(i, j) for i in range(m) for j in range(n) if C[i, j] > 0]
    zero_costs = [(i, j) for i in range(m) for j in range(n) if C[i, j] == 0]

    positive_costs.sort(key=lambda pos: (C[pos[0], pos[1]], pos[0], pos[1]))
    zero_costs.sort(key=lambda pos: (pos[0], pos[1]))

    for cost_list in [positive_costs, zero_costs]:
        for i, j in cost_list:
            if a[i] == 0 or b[j] == 0:
                continue
            value = min(a[i], b[j])
            x[i, j] = value
            a[i] -= value
            b[j] -= value

    return x

def print_iteration_table(title, C, X):
    headers = [f"B{j+1}" for j in range(X.shape[1])] + ["Σ"]
    table = []
    for i in range(X.shape[0]):
        row = []
        for j in range(X.shape[1]):
            value = int(X[i][j]) if X[i][j] > 0 else "-"
            cost = int(C[i][j])
            row.append(f"{cost}[{value}]" if value != "-" else str(cost))
        row_sum = int(np.sum(X[i, :]))
        table.append([f"A{i+1}"] + row + [row_sum])
    footer = ["Σ"] + [int(np.sum(X[:, j])) for j in range(X.shape[1])] + [""]
    print(f"\n{title}")
    print(tabulate(table + [footer], headers=[""] + headers, tablefmt="fancy_grid"))

def method_potentials(supply, demand, costs, initial_plan):
    C = np.copy(costs)
    X = initial_plan.astype(float)
    n, m = C.shape

    iteration = 0
    while True:
        iteration += 1
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

        print_iteration_table(f"\n--- Ітерація {iteration} ---", C, X)

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

def print_cost_expression(costs, plan):
    expr = []
    total = 0
    for i in range(costs.shape[0]):
        for j in range(costs.shape[1]):
            if plan[i, j] > 0:
                c, x = int(costs[i, j]), int(plan[i, j])
                expr.append(f"{c}*{x}")
                total += c * x
    print(">>> Обчислення загальної вартості перевезення")
    print("Z = " + " + ".join(expr) + f" = {total}\n")

if __name__ == '__main__':
    supply = np.array([130, 90, 40])
    demand = np.array([110, 30, 50, 80, 90])
    costs = np.array([[4, 5, 6, 8, 10], [10, 3, 2, 3, 15], [4, 10, 5, 1, 16]])

    supply, demand, costs, is_open, added_type = check_problem_type(supply, demand, costs)

    print_transport_table(costs, supply, demand, "Вартість перевезень")

    print("\n--- Вибір методу ---")
    print("[1] Північно-західний кут")
    print("[2] Найменша вартість")
    method = input("Введіть номер методу (1 або 2): ")

    if method == '1':
        initial_plan = northwest_corner_method(supply, demand)
        method_name = "Північно-західний кут"
    elif method == '2':
        initial_plan = least_cost_method(costs, supply, demand)
        method_name = "Найменша вартість"
    else:
        print("Некоректний вибір методу.")
        exit()

    print_transport_table(costs, supply, demand, f"Початковий план — {method_name}", initial_plan)
    print_cost_expression(costs, initial_plan)

    plan, total_cost = method_potentials(supply, demand, costs, initial_plan)

    if plan is not None:
        print_transport_table(costs[:len(supply), :len(demand)], supply, demand, "Оптимізований план", np.round(plan).astype(int))
        print_cost_expression(costs[:len(supply), :len(demand)], np.round(plan).astype(int))

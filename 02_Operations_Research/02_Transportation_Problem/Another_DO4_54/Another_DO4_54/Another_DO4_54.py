import numpy as np
from collections import Counter
from tabulate import tabulate

def nw_corner_rule(supply, demand):
    s, d = supply.copy(), demand.copy()
    rows, cols = len(s), len(d)
    plan = np.full((rows, cols), np.nan)
    row_idx = col_idx = 0
    while row_idx < rows and col_idx < cols:
        allocation = min(s[row_idx], d[col_idx])
        plan[row_idx, col_idx] = allocation
        s[row_idx] -= allocation
        d[col_idx] -= allocation
        if s[row_idx] == 0 and d[col_idx] == 0:
            if col_idx + 1 < cols:
                col_idx += 1
            elif row_idx + 1 < rows:
                row_idx += 1
            else:
                break
        elif s[row_idx] == 0:
            row_idx += 1
        elif d[col_idx] == 0:
            col_idx += 1
    return np.nan_to_num(plan, nan=0).astype(int)

def cheapest_value(cost_matrix, supply, demand):
    s, d = supply.copy(), demand.copy()
    rows, cols = len(s), len(d)
    plan = np.zeros((rows, cols), dtype=int)
    positive_cells = [(i, j) for i in range(rows) for j in range(cols) if cost_matrix[i, j] > 0]
    zero_cells = [(i, j) for i in range(rows) for j in range(cols) if cost_matrix[i, j] == 0]
    
    positive_cells.sort(key=lambda pos: (cost_matrix[pos], pos[0], pos[1]))
    zero_cells.sort(key=lambda pos: (pos[0], pos[1]))

    for cell_group in [positive_cells, zero_cells]:
        for r, c in cell_group:
            if s[r] == 0 or d[c] == 0:
                continue
            quantity = min(s[r], d[c])
            plan[r, c] = quantity
            s[r] -= quantity
            d[c] -= quantity
    return plan

def optimizer(supply, demand, costs, plan):
    C, X = np.copy(costs), plan.astype(float)
    n, m = C.shape
    while True:
        u, v= np.array([np.nan] * n), np.array([np.nan] * m)
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

def compute_total_cost(c, x):
    total = sum(c[i,j]*x[i,j] for i in range(c.shape[0]) for j in range(c.shape[1]))
    print("\nОптимум цільвої F(x)" + f" = {int(total)}\n")

def display_table(C, supply, demand, plan=None, heading=""):
    headers = [f"B{j+1}" for j in range(len(demand))] + ["Available"]
    matrix = []

    for i, s in enumerate(supply):
        row = [f"A{i+1}"]
        for j, d in enumerate(demand):
            if plan is not None and plan[i, j] > 0:
                cell = f"{int(C[i,j])}[{int(plan[i,j])}]"
            else:
                cell = "-" if C[i,j] == 0 else str(int(C[i,j]))
            row.append(cell)
        row.append(int(s))
        matrix.append(row)

    footer = ["Demand"] + [int(x) for x in demand] + [""]
    print(f"\n{heading}")
    print(tabulate(matrix + [footer], headers=[""] + headers, tablefmt="presto"))

def balance_problem(supply, demand, costs):
    s_sum, d_sum = sum(supply), sum(demand)
    print(f"→ Перевірка балансу: Запаси = {s_sum}, Потреби = {d_sum}")

    if s_sum == d_sum:
        print("✓ Задача збалансована.")
        return supply, demand, costs, False, None

    print("⚠ Задача незбалансована.")
    difference = abs(s_sum - d_sum)

    if s_sum < d_sum:
        print(f"→ Додаємо фіктивного постачальника із запасом {difference}")
        supply = np.append(supply, difference)
        costs = np.vstack([costs, np.zeros((1, costs.shape[1]))])
        return supply, demand, costs, True, 'supplier'
    else:
        print(f"→ Додаємо фіктивного споживача із потребою {difference}")
        demand = np.append(demand, difference)
        costs = np.hstack([costs, np.zeros((costs.shape[0], 1))])
        return supply, demand, costs, True, 'consumer'

def start():
    print("\nОберіть метод побудови початкового плану:")
    print("[1] Правило північно-західного кута (NWC)")
    print("[2] Метод найменших витрат (LC)")

    choice = input("\nВвід: ").strip()
    if choice == "1":
        method_choice = "NWC"
    elif choice == "2":
        method_choice = "LC"
    else:
        print("Некоректний вибір. Завершення програми.")
        return

    supply = np.array([50, 20, 30, 20])
    demand = np.array([50, 30, 35, 15])
    costs = np.array([ [1, 3, 3, 4], [5, 2, 7, 5], [6, 4, 8, 2],[7, 1, 5, 7]])
    print("\n--- Перевірка збалансованості ---")
    supply, demand, costs, *_ = balance_problem(supply, demand, costs)

    if method_choice == "NWC":
        initial_plan, label = nw_corner_rule(supply, demand), "Правило північно-західного кута"
    else:
        initial_plan, label = cheapest_value(costs, supply, demand), "Метод найменших витрат"

    print(f"\nПочатковий план ({label}):")
    display_table(costs, supply, demand, plan=initial_plan)
    compute_total_cost(costs, initial_plan)
    final_plan, final_cost = optimizer(supply, demand, costs, initial_plan)

    print("\nОптимізований план (методом потенціалів):")
    display_table(costs[:len(supply), :len(demand)], supply, demand,
                  plan=np.round(final_plan).astype(int))
    compute_total_cost(costs[:len(supply), :len(demand)], np.round(final_plan).astype(int))

if __name__ == "__main__":
    start()
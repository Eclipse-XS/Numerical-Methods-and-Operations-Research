import numpy as np
from collections import Counter
from tabulate import tabulate

def fancy_title(text, char="=", width=50):
    line = char * width
    print(f"\n{line}\n{text.center(width)}\n{line}")

def print_transport_table(c, a, b, title="Транспортна таблиця", plan=None):
    hdr = [f"B{j+1}" for j in range(len(b))] + ["aᵢ"]
    tbl = [[f"A{i+1}"] + [f"{int(c[i,j])}[{int(plan[i,j])}]" if plan is not None and plan[i,j] > 0 else "-" if c[i,j] == 0 else str(int(c[i,j])) for j in range(len(b))] + [int(a[i])] for i in range(len(a))]
    print(f"\n{title}")
    print(tabulate(tbl + [["bⱼ"] + list(map(int, b)) + [""]], headers=[""] + hdr, tablefmt="github"))

def check_problem_type(a, b, c):
    fancy_title("ТРАНСПОРТНА ЗАДАЧА", "=")
    print("Тип задачі:")
    
    sa, sb = sum(a), sum(b)
    
    if sa == sb:
        print("- Задача є збалансованою.")
        return a, b, c, False, None

    diff = abs(sa - sb)

    print("- Задача є незбалансованою.")
    
    if sa < sb:
        print(f"- ∑a = {sa} < ∑b = {sb} → додаємо фіктивного **постачальника** з запасом `a = {diff}`.")
        a = np.append(a, diff)
        c = np.vstack([c, np.zeros((1, c.shape[1]))])
        return a, b, c, True, 'supplier'
    else:
        print(f"- ∑a = {sa} > ∑b = {sb} → додаємо фіктивного **споживача** з потребою `b = {diff}`.")
        b = np.append(b, diff)
        c = np.hstack([c, np.zeros((c.shape[0], 1))])
        return a, b, c, True, 'consumer'


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

def method_potentials(supply, demand, costs, initial_plan):
    C = np.copy(costs)
    X = initial_plan.astype(float)
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

def print_cost_expression(c, x):
    terms = [f"{int(c[i,j])}*{int(x[i,j])}" for i in range(c.shape[0]) for j in range(c.shape[1]) if x[i,j] > 0]
    total = sum(c[i,j]*x[i,j] for i in range(c.shape[0]) for j in range(c.shape[1]))
    print("\nРозрахунок загальної вартості:")
    print("Z = " + " + ".join(terms) + f" = {int(total)} грн\n")

if __name__ == '__main__':
    METHOD = "LC"  # ← змінюй тут на "NWC" або "LC"

    supply = np.array([20, 20, 40, 45])
    demand = np.array([25, 30, 40, 30])
    costs = np.array([[1, 3, 3, 8], [8, 6, 2, 6], [7, 7, 3, 8], [5, 2, 4, 5]])

    supply, demand, costs, *_ = check_problem_type(supply, demand, costs)

    if METHOD == "NWC":
        plan0, name = northwest_corner_method(supply, demand), "Північно-західний кут"
    elif METHOD == "LC":
        plan0, name = least_cost_method(costs, supply, demand), "Метод найменшої вартості"
    else:
        print("Невідомий метод."); exit()

    fancy_title(f"Початковий план — {name}", "-")
    print_transport_table(costs, supply, demand, plan=plan0)
    print_cost_expression(costs, plan0)

    final_plan, total_cost = method_potentials(supply, demand, costs, plan0)
    fancy_title("Оптимізований план", "=")
    print_transport_table(costs[:len(supply), :len(demand)], supply, demand, plan=np.round(final_plan).astype(int))
    print_cost_expression(costs[:len(supply), :len(demand)], np.round(final_plan).astype(int))
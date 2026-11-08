import numpy as np
from collections import Counter
from tabulate import tabulate

def print_cost_expression(c, x):
    terms = [f"{int(c[i,j])}*{int(x[i,j])}" for i in range(c.shape[0]) for j in range(c.shape[1]) if x[i,j] > 0]
    total = sum(c[i,j]*x[i,j] for i in range(c.shape[0]) for j in range(c.shape[1]))
    print("\nМінімум функції:")
    print("F(x)" +  f" = {int(total)} грн\n")

def fancy_title(text, char="=", width=50):
    line = char * width
    print(f"\n{line}\n{text.center(width)}\n{line}")

def print_transport_table(c, a, b, title="Транспортна таблиця", plan=None, dummy_info=None):
    dummy = dummy_info[0] if dummy_info else False
    dummy_type = dummy_info[1] if dummy_info else None

    hdr = [f"B{j+1}" for j in range(len(b))]
    if dummy and dummy_type == 'consumer':
        hdr[-1] += " (фікт.)"
    hdr += ["aᵢ"]

    tbl = []
    for i in range(len(a)):
        row_label = f"A{i+1}"
        if dummy and dummy_type == 'supplier' and i == len(a) - 1:
            row_label += " (фікт.)"

        row = [row_label]
        for j in range(len(b)):
            cij = int(c[i, j])
            if plan is not None and plan[i, j] > 0:
                cell = f"{cij}[{int(plan[i, j])}]"
            elif cij == 0 and dummy and ((dummy_type == 'supplier' and i == len(a) - 1) or (dummy_type == 'consumer' and j == len(b) - 1)):
                cell = "-"
            else:
                cell = str(cij)
            row.append(cell)
        row.append(int(a[i]))
        tbl.append(row)

    footer = ["bⱼ"] + list(map(int, b)) + [""]

    print(f"\n{title}")
    print(tabulate(tbl + [footer], headers=[""] + hdr, tablefmt="pretty"))


def check_problem_type(a, b, c):
    sa, sb = sum(a), sum(b)

    print(f"∑a (запаси)  = {sa}")
    print(f"∑b (потреби) = {sb}\n")

    if sa == sb:
        print("Задача збалансована.")
        return a, b, c, False, None

    diff = abs(sa - sb)
    print("Задача є незбалансованою.")

    if sa < sb:
        print(f"Вводимо фіктивного постачальника із запасом a = {diff}.")
        a = np.append(a, diff)
        c = np.vstack([c, np.zeros((1, c.shape[1]))])
        return a, b, c, True, 'supplier'
    else:
        print(f"Вводимо фіктивного споживача із потребою b = {diff}.")
        b = np.append(b, diff)
        c = np.hstack([c, np.zeros((c.shape[0], 1))])
        return a, b, c, True, 'consumer'

def nw_corner(a, b):
    a, b = a.copy(), b.copy()
    m, n = len(a), len(b)
    x = np.zeros((m, n), dtype=int)
    i = j = 0
    while i < m and j < n:
        v = min(a[i], b[j])
        x[i, j] = v
        a[i] -= v
        b[j] -= v
        if a[i] == 0 and b[j] == 0:
            if j + 1 < n:
                x[i, j+1] = 0
                j += 1
            elif i + 1 < m:
                x[i+1, j] = 0
                i += 1
            else:
                break
        elif a[i] == 0:
            i += 1
        elif b[j] == 0:
            j += 1
    return x

def least_cost(C, a, b):
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

def potentials(supply, demand, costs, initial_plan):
    C, X = np.copy(costs), initial_plan.astype(float)
    n, m = C.shape
    while True:
        u, v = np.full(n, np.nan), np.full(m, np.nan)
        S = np.zeros((n, m))
        _x, _y = np.where(X > 0)
        nonzero = list(zip(_x, _y))
        if len(nonzero) < n + m - 1:
            for i in range(n):
                for j in range(m):
                    if X[i, j] == 0:
                        X[i, j] = 1e-5
                        nonzero.append((i, j))
                        if len(nonzero) >= n + m - 1: break
                if len(nonzero) >= n + m - 1: break
        u[nonzero[0][0]] = 0
        changed = True
        while changed:
            changed = False
            for i, j in nonzero:
                if np.isnan(u[i]) and not np.isnan(v[j]):
                    u[i], changed = C[i, j] - v[j], True
                elif not np.isnan(u[i]) and np.isnan(v[j]):
                    v[j], changed = C[i, j] - u[i], True

        S = C - u[:, None] - v
        if np.all(S >= 0):
            break
        i, j = np.unravel_index(np.argmin(S), S.shape)
        start = (i, j)
        T = np.copy(X)
        T[start] = 1
        for _ in range(1000):
            _xs, _ys = np.nonzero(T)
            xcount, ycount = Counter(_xs), Counter(_ys)
            for x in [k for k, v in xcount.items() if v <= 1]:
                T[x, :] = 0
            for y in [k for k, v in ycount.items() if v <= 1]:
                T[:, y] = 0
            if all(v > 1 for v in xcount.values()) and all(v > 1 for v in ycount.values()):
                break
        else:
            return X, np.sum(X * C)
        dist = lambda a, b: abs(a[0] - b[0]) + abs(a[1] - b[1])
        fringe = set(map(tuple, np.argwhere(T > 0)))
        path = [start]
        while len(path) < len(fringe):
            last = path[-1]
            fringe.discard(last)
            candidates = [p for p in fringe if p[0] == last[0] or p[1] == last[1]]
            if not candidates: return X, np.sum(X * C)
            path.append(min(candidates, key=lambda x: dist(last, x)))
        neg, pos = path[1::2], path[::2]
        q = min(X[i, j] for i, j in neg)
        for i, j in neg: X[i, j] -= q
        for i, j in pos: X[i, j] += q
    return X, np.sum(X * C)

METHOD = "NWC"

ai = np.array([30, 70, 50])
bj = np.array([30, 40, 20, 60])
cij = np.array([
    [2, 7, 3, 6],
    [9, 4, 5, 7],
    [5, 7, 6, 2]
])

ai, bj, cij, is_dummy, dummy_type = check_problem_type(ai, bj, cij)

if METHOD == "NWC":
    plan0 = nw_corner(ai, bj)
    method_name = "Північно-західний кут"
elif METHOD == "LC":
    plan0 = least_cost(cij, ai, bj)
    method_name = "Метод найменшої вартості"
else:
    raise ValueError(f"Невідомий метод: {METHOD}")

print_transport_table(cij, ai, bj,f"{method_name}", plan=plan0, dummy_info=(is_dummy, dummy_type))
print_cost_expression(cij, plan0)

final_plan, total_cost = potentials(ai, bj, cij, plan0)

print_transport_table(cij[:len(ai), :len(bj)], ai, bj, "Метод потенціалів", plan=np.round(final_plan).astype(int), dummy_info=(is_dummy, dummy_type))
print_cost_expression(cij[:len(ai), :len(bj)], np.round(final_plan).astype(int))
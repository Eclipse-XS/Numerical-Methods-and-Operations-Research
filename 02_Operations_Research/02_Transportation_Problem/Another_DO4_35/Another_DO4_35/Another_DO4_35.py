import os
import numpy as np
from collections import Counter

def NWC(supply, demand):
    supply = supply.copy()
    demand = demand.copy()
    num_rows, num_cols = len(supply), len(demand)
    allocation = np.full((num_rows, num_cols), np.nan)
    i = j = 0

    while i < num_rows and j < num_cols:
        qty = min(supply[i], demand[j])
        allocation[i, j] = qty
        supply[i] -= qty
        demand[j] -= qty

        if supply[i] == 0 and demand[j] == 0:
            if j + 1 < num_cols:
                j += 1
            elif i + 1 < num_rows:
                i += 1
            else:
                break
        elif supply[i] == 0:
            i += 1
        elif demand[j] == 0:
            j += 1

    return np.nan_to_num(allocation, nan=0).astype(int)

def LC(costs, supply, demand):
    supply = supply.copy()
    demand = demand.copy()
    num_rows, num_cols = len(supply), len(demand)
    allocation = np.zeros((num_rows, num_cols), dtype=int)

    pos_cells = [(i, j) for i in range(num_rows) for j in range(num_cols) if costs[i, j] > 0]
    zero_cells = [(i, j) for i in range(num_rows) for j in range(num_cols) if costs[i, j] == 0]

    pos_cells.sort(key=lambda cell: (costs[cell[0], cell[1]], cell[0], cell[1]))
    zero_cells.sort(key=lambda cell: (cell[0], cell[1]))

    for cell_list in [pos_cells, zero_cells]:
        for i, j in cell_list:
            if supply[i] == 0 or demand[j] == 0:
                continue
            qty = min(supply[i], demand[j])
            allocation[i, j] = qty
            supply[i] -= qty
            demand[j] -= qty

    return allocation

def MP(supply, demand, costs, plan):
    cost_matrix = np.copy(costs)
    allocation = plan.astype(float)
    num_rows, num_cols = cost_matrix.shape

    while True:
        u = np.full(num_rows, np.nan)
        v = np.full(num_cols, np.nan)
        reduced_costs = np.zeros((num_rows, num_cols))

        filled_cells = list(zip(*np.where(allocation > 0)))

        if len(filled_cells) < num_rows + num_cols - 1:
            for i in range(num_rows):
                for j in range(num_cols):
                    if allocation[i, j] == 0:
                        allocation[i, j] = 0.00001
                        filled_cells.append((i, j))
                        if len(filled_cells) >= num_rows + num_cols - 1:
                            break
                if len(filled_cells) >= num_rows + num_cols - 1:
                    break

        u[filled_cells[0][0]] = 0
        updated = True
        while updated:
            updated = False
            for i, j in filled_cells:
                if np.isnan(u[i]) and not np.isnan(v[j]):
                    u[i] = cost_matrix[i, j] - v[j]
                    updated = True
                elif not np.isnan(u[i]) and np.isnan(v[j]):
                    v[j] = cost_matrix[i, j] - u[i]
                    updated = True

        for i in range(num_rows):
            for j in range(num_cols):
                reduced_costs[i, j] = cost_matrix[i, j] - u[i] - v[j]

        if np.all(reduced_costs >= 0):
            break

        pivot_row, pivot_col = np.unravel_index(np.argmin(reduced_costs), reduced_costs.shape)
        start = (pivot_row, pivot_col)
        T = np.copy(allocation)
        T[start] = 1

        for _ in range(1000):
            row_counts = Counter(np.nonzero(T)[0])
            col_counts = Counter(np.nonzero(T)[1])
            for i, count in row_counts.items():
                if count <= 1:
                    T[i, :] = 0
            for j, count in col_counts.items():
                if count <= 1:
                    T[:, j] = 0
            if all(c > 1 for c in row_counts.values()) and all(c > 1 for c in col_counts.values()):
                break

        def manhattan(p1, p2): return abs(p1[0] - p2[0]) + abs(p1[1] - p2[1])
        fringe = set(map(tuple, np.argwhere(T > 0)))
        path = [start]
        while len(path) < len(fringe):
            last = path[-1]
            if last in fringe:
                fringe.remove(last)
            candidates = [cell for cell in fringe if cell[0] == last[0] or cell[1] == last[1]]
            if not candidates:
                return allocation, np.sum(allocation * cost_matrix)
            next_step = min(candidates, key=lambda cell: manhattan(last, cell))
            path.append(next_step)

        negative_cells = path[1::2]
        positive_cells = path[::2]
        min_value = min(allocation[i, j] for i, j in negative_cells)

        for i, j in negative_cells:
            allocation[i, j] -= min_value
        for i, j in positive_cells:
            allocation[i, j] += min_value

    return allocation, np.sum(allocation * cost_matrix)


def clear_screen():
    os.system("cls" if os.name == "nt" else "clear")

def fancy_title(text, char="=", width=50):
    line = char * width
    print(f"\n{line}\n{text.center(width)}\n{line}")

def print_custom_table(c, a, b, title="Transport Table", plan=None):
    m, n = len(a), len(b)
    col_headers = [f"B{j+1}" for j in range(n)]
    row_headers = [f"A{i+1}" for i in range(m)]

    print(f"\n{title}")
    print("      " + "  ".join(f"{h:^8}" for h in col_headers) + "    Supply")

    for i in range(m):
        row = [f"{row_headers[i]}"]
        for j in range(n):
            if plan is not None and plan[i][j] > 0:
                val = f"{int(c[i][j])}[{int(plan[i][j])}]"
            elif c[i][j] == 0:
                val = "-"
            else:
                val = f"{int(c[i][j])}"
            row.append(f"{val:^8}")
        row.append(f"{int(a[i]):>6}")
        print("  ".join(row))

    demand_row = ["Demand"] + [f"{int(b[j]):^8}" for j in range(n)] + [""]
    print("  ".join(demand_row))

def check_problem_type(a, b, c):
    fancy_title("TRANSPORTATION PROBLEM", "=")
    print("Problem type:")

    sa, sb = sum(a), sum(b)

    if sa == sb:
        print("- Problem is balanced.")
        return a, b, c, False, None

    diff = abs(sa - sb)
    print("- Problem is unbalanced.")

    if sa < sb:
        print(f"- Total supply = {sa} < Total demand = {sb} → Adding dummy supplier with supply = {diff}.")
        a = np.append(a, diff)
        c = np.vstack([c, np.zeros((1, c.shape[1]))])
        return a, b, c, True, 'supplier'
    else:
        print(f"- Total supply = {sa} > Total demand = {sb} → Adding dummy consumer with demand = {diff}.")
        b = np.append(b, diff)
        c = np.hstack([c, np.zeros((c.shape[0], 1))])
        return a, b, c, True, 'consumer'

def print_cost_expression(c, x):
    total = int(np.sum(c * x))
    print("\nCost Calculation:")
    print(f"F(x) --> min = {total} $\n")

def run_transportation_problem():
    print("Select initial basic feasible solution method:")
    print("1. Northwest Corner Method")
    print("2. Least Cost Method")
    method_choice = input("Enter 1 or 2: ").strip()

    clear_screen()

    method_name = "Northwest Corner"
    supply = np.array([50, 40, 20])
    demand = np.array([30, 25, 35, 20])
    costs = np.array([
        [3, 2, 4, 1],
        [2, 3, 1, 5],
        [3, 2, 4, 4]
    ])

    supply, demand, costs, *_ = check_problem_type(supply, demand, costs)

    if method_choice == "1":
        initial_plan = NWC(supply, demand)
    elif method_choice == "2":
        initial_plan = LC(costs, supply, demand)
        method_name = "Least Cost"
    else:
        print("Unknown choice. Exiting.")
        return

    fancy_title(f"Initial Plan — {method_name}", "-")
    print_custom_table(costs, supply, demand, plan=initial_plan)
    print_cost_expression(costs, initial_plan)

    optimized_plan, total_cost = MP(supply, demand, costs, initial_plan)

    fancy_title("Optimized Plan", "=")
    print_custom_table(costs[:len(supply), :len(demand)], supply, demand, plan=np.round(optimized_plan).astype(int))
    print_cost_expression(costs[:len(supply), :len(demand)], np.round(optimized_plan).astype(int))

run_transportation_problem()
import numpy as np
from tabulate import tabulate

profits = np.array([
    [0,   0,   0,   0],
    [50,  51,  52,  50],
    [60,  59,  59,  55],
    [67,  65,  65,  63],
    [73,  72,  73,  72]
])

num_enterprises = 4
total_investment = 40
step = 10
steps = total_investment // step + 1

F = np.zeros((num_enterprises + 2, steps), dtype=int)
decision = np.zeros((num_enterprises + 1, steps), dtype=int)

# Динамічне програмування
for i in range(num_enterprises, 0, -1):
    for x in range(steps):
        max_val = -1
        best_k = 0
        for k in range(x + 1):
            current_profit = profits[k][i - 1] + F[i + 1][x - k]
            if current_profit > max_val:
                max_val = current_profit
                best_k = k
        F[i][x] = max_val
        decision[i][x] = best_k

# Вивід таблички з tabulate
def print_stage_tabulate(i):
    xi = f"x{i} (тис. грн)"
    ki_labels = [f"k{i} = {k * step}" for k in range(steps)]
    fi = f"f{i}(x{i})"
    headers = [xi] + ki_labels + [fi, f"k{i}*"]
    
    table = []
    for x in range(steps):
        row = [x * step]
        max_val = -1
        best_k = 0
        for k in range(steps):
            if k <= x:
                a = profits[k][i - 1]
                b = F[i + 1][x - k]
                val = a + b
                expr = f"{a} + {b} = {val}" if i < 4 else str(a)
                row.append(expr)
                if val > max_val:
                    max_val = val
                    best_k = k * step
            else:
                row.append("")
        row.append(max_val)
        row.append(best_k)
        table.append(row)

    print(f"\n=== ЕТАП {i} ===")
    print(tabulate(table, headers=headers, tablefmt="fancy_grid"))

# Вивід етапів
for i in range(num_enterprises, 0, -1):
    print_stage_tabulate(i)

# Відновлення оптимального плану
remaining = total_investment // step
allocation = [0] * num_enterprises
for i in range(1, num_enterprises + 1):
    alloc_k = decision[i][remaining]
    allocation[i - 1] = alloc_k * step
    remaining -= alloc_k

# Вивід результату
print("\n=== ОПТИМАЛЬНИЙ ПЛАН ІНВЕСТИЦІЙ ===")
total_profit = 0
summary = []
for i in range(num_enterprises):
    inv = allocation[i]
    k_index = inv // step
    profit = profits[k_index][i]  # КОРЕКТНА індексація
    total_profit += profit
    summary.append([f"Підприємство {i + 1}", f"{inv} тис. грн", f"{profit} тис. грн"])

print(tabulate(summary, headers=["Підприємство", "Інвестиція", "Прибуток"], tablefmt="grid"))
print(f"\nСумарний прибуток: {total_profit} тис. грн")

import numpy as np
from tabulate import tabulate
from colorama import init, Fore, Style

init(autoreset=True)

profits = np.array([
    [0,   0,   0,   0],
    [115, 135, 145, 105],
    [315, 285, 295, 215],
    [445, 405, 405, 465],
    [565, 545, 575, 555],
    [695, 715, 715, 695]
])

num_enterprises = 4
total_investment = 1000
step = 200
steps = total_investment // step + 1

F = np.zeros((num_enterprises + 2, steps), dtype=int)
decision = np.zeros((num_enterprises + 1, steps), dtype=int)

output = []

for i in range(num_enterprises, 0, -1):
    results_table = []
    headers = [f"{k*step} грн" for k in range(steps)] + ["🔝 Макс. прибуток", "💸 Оптим. вкладення"]
    for x in range(steps):
        row = []
        max_val = -1
        best_k = 0
        for k in range(steps):
            if k <= x:
                val = profits[k][i - 1] + F[i + 1][x - k]
                row.append(val)
                if val > max_val:
                    max_val = val
                    best_k = k * step
            else:
                row.append("")
        row.append(max_val)
        row.append(best_k)
        results_table.append(row)

    table_str = tabulate(
        results_table,
        headers=headers,
        showindex=[f"{x*step} грн" for x in range(steps)],
        tablefmt="github"
    )

    stage_title = f"{Fore.CYAN}{Style.BRIGHT}📊 ЕТАП {i} — Розподіл для підприємства {i}{Style.RESET_ALL}"
    output.append(stage_title)
    output.append(table_str)

    for x in range(steps):
        F[i][x] = results_table[x][-2]
        decision[i][x] = results_table[x][-1] // step

remaining = total_investment // step
allocation = [0] * num_enterprises

for i in range(1, num_enterprises + 1):
    alloc_k = decision[i][remaining]
    allocation[i - 1] = alloc_k * step
    remaining -= alloc_k

output.append(f"\n{Fore.GREEN}{Style.BRIGHT}✅ ОПТИМАЛЬНИЙ РОЗПОДІЛ ІНВЕСТИЦІЙ:{Style.RESET_ALL}")

table = []
total_profit = 0

for i in range(num_enterprises):
    investment = allocation[i]
    profit = profits[investment // step][i]
    table.append([f"🏢 Підприємство {i+1}", f"{investment} тис. грн", f"{profit} тис. грн"])
    total_profit += profit

output.append(tabulate(table, headers=["Підприємство", "Інвестиція", "Прибуток"], tablefmt="github"))
output.append(tabulate([["💰 Сумарний прибуток", f"{total_profit} тис. грн"]], tablefmt="github"))

full_output = "\n\n".join(output)
print(full_output)

import numpy as np

profits = np.array([
    [0, 0, 0, 0],
    [30, 31, 32, 30],
    [40, 39, 39, 35],
    [47, 45, 45, 43],
    [53, 52, 53, 52],
    [59, 60, 60, 59]
])

num_enterprises = 4
total_investment = 50
step = 10
steps = total_investment // step + 1

F = np.zeros((num_enterprises + 2, steps), dtype=int)
decision = np.zeros((num_enterprises + 1, steps), dtype=int)

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

def print_stage(i):
    print(f"\n=== ЕТАП {i} (ПІДПРИЄМСТВО {i}) ===")
    header = f"{'x'+str(i)+' (тис. грн)':<18}" + "".join([f"{k*step:>8}" for k in range(steps)])
    header += f"{'|':>3}{'f_'+str(i)+'(x)':>10}{'k*':>8}"
    print(header)
    print("-" * len(header))
    
    for x in range(steps):
        row = f"{x * step:<18}"
        max_val = -1
        best_k = 0
        for k in range(steps):
            if k <= x:
                val = profits[k][i - 1] + F[i + 1][x - k]
                row += f"{val:>8}"
                if val > max_val:
                    max_val = val
                    best_k = k * step
            else:
                row += f"{'':>8}"
        row += f"{'|':>3}{max_val:>10}{best_k:>8}"
        print(row)


for i in range(num_enterprises, 0, -1):
    print_stage(i)

remaining = total_investment // step
allocation = [0] * num_enterprises

for i in range(1, num_enterprises + 1):
    alloc_k = decision[i][remaining]
    allocation[i - 1] = alloc_k * step
    remaining -= alloc_k

print("\n========== ОПТИМАЛЬНИЙ РОЗВ’ЯЗОК ==========")
total_profit = 0
for i in range(num_enterprises):
    p = profits[allocation[i] // step][i]
    print(f"  Підприємство {i+1}: інвестиція = {allocation[i]} тис. грн, прибуток = {p} тис. грн")
    total_profit += p
print(f"\n  Сумарний прибуток: {total_profit} тис. грн")

print("=" * 40)
print("Група: ОІ-23")
print("Виконав: Гончаров Вадим Олегович")
print("Номер ЛР: №7")
print("Варіант: №4")
print("Номер студентського квитка: №75")
print("=" * 40)
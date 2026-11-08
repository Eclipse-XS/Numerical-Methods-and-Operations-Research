import numpy as np
from tabulate import tabulate

def evaluate_optimal_values(total_firms, num_levels, revenue_matrix, value_matrix, decision_matrix):
    for firm in range(total_firms, 0, -1):
        for budget_level in range(num_levels):
            max_result = -1
            selected_step = 0
            for step_idx in range(num_levels):
                if step_idx <= budget_level:
                    candidate_value = revenue_matrix[step_idx][firm - 1] + value_matrix[firm + 1][budget_level - step_idx]
                    if candidate_value > max_result:
                        max_result = candidate_value
                        selected_step = step_idx
            value_matrix[firm][budget_level] = max_result
            decision_matrix[firm][budget_level] = selected_step

def show_stage_table(firm_index, total_firms, num_levels, step_unit, revenue_matrix, value_matrix):
    header_row = [f"x{firm_index} (тис. грн)"] + [f"k = {k * step_unit}" for k in range(num_levels)] + [f"f{firm_index}(x)", f"k{firm_index}*"]
    output_table = []

    for level in range(num_levels):
        row_data = [level * step_unit]
        max_score = -1
        best_step_val = 0
        for s in range(num_levels):
            if s <= level:
                result = revenue_matrix[s][firm_index - 1] + value_matrix[firm_index + 1][level - s]
                row_data.append(result)
                if result > max_score:
                    max_score = result
                    best_step_val = s * step_unit
            else:
                row_data.append("")
        row_data.append(max_score)
        row_data.append(best_step_val)
        output_table.append(row_data)

    print(f"\n=== СТАДІЯ {firm_index} (ФІРМА {firm_index}) ===")
    print(tabulate(output_table, headers=header_row, tablefmt="presto"))


def restore_optimal_plan(total_firms, step_unit, total_budget, decision_matrix):
    allocation = []
    remaining = total_budget // step_unit
    for i in range(1, total_firms + 1):
        step = decision_matrix[i][remaining]
        allocation.append(step * step_unit)
        remaining -= step
    return allocation

def present_final_result(allocation, revenue_matrix, step_unit):
    total_income = 0
    print("\n=== ОПТИМАЛЬНИЙ РОЗПОДІЛ ІНВЕСТИЦІЙ ===\n")
    for i, invested in enumerate(allocation):
        income = revenue_matrix[invested // step_unit][i]
        total_income += income
        print(f"- Підприємство {i + 1} — {invested} тис. грн → прибуток {income} тис. грн")
    print(f"\nМаксимальний сумарний прибуток становить {total_income} тис. грн.")

def main():
    revenue_matrix = np.array([
        [0,    0,    0,    0],
        [120,  140,  150,  110],
        [320,  290,  300,  220],
        [450,  410,  410,  470],
        [570,  550,  580,  560],
        [700,  720,  720,  700]
    ])

    total_firms = 4
    total_budget = 1000
    step_unit = 200
    num_levels = total_budget // step_unit + 1

    value_matrix = np.zeros((total_firms + 2, num_levels), dtype=int)
    decision_matrix = np.zeros((total_firms + 1, num_levels), dtype=int)

    evaluate_optimal_values(total_firms, num_levels, revenue_matrix, value_matrix, decision_matrix)

    for firm_index in range(total_firms, 0, -1):
        show_stage_table(firm_index, total_firms, num_levels, step_unit, revenue_matrix, value_matrix)

    allocation = restore_optimal_plan(total_firms, step_unit, total_budget, decision_matrix)
    present_final_result(allocation, revenue_matrix, step_unit)

if __name__ == "__main__":
    main()
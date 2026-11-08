import numpy as np
import matplotlib.pyplot as plt

h = 0.01 
x0, x_end = 0, 2 
y0 = 0.5 

def f(x, y):
    return y * x**2 + 3 * np.sin(y)

print("Виберіть метод обчислення:")
print("1 - Метод Ейлера")
print("2 - Метод Рунге-Кутта")
choice = input("Введіть 1 або 2: ")

x_values = np.arange(x0, x_end + h, h)
y_values = [y0]

if choice == '1':
    f_values = []
    hf_values = []

    for i in range(1, len(x_values)):
        x_prev, y_prev = x_values[i - 1], y_values[-1]
        f_val = f(x_prev, y_prev)
        # Метод Ейлера
        y_next = y_prev + h * f_val

        y_values.append(y_next)
        f_values.append(f_val)
        hf_values.append(h * f_val)

    print(f"\nМетод Ейлера:\n{'i':<5}{'x_i':<10}{'y_i':<15}{'f(x_i, y_i)':<20}{'hf(x_i, y_i)':<20}")
    print("-" * 70)
    for i in range(len(x_values)):
        x = x_values[i]
        y = y_values[i]
        f_val = f_values[i - 1] if i >= 0 else 0.0
        hf_val = hf_values[i - 1] if i >= 0 else 0.0
        
        print(f"{i:<5}{x:<10.2f}{y:<15.5f}{f_val:<20.10f}{hf_val:<20.10f}")

    plt.figure(figsize=(10, 6))
    plt.plot(x_values, y_values, label="Наближений вигляд інтегральної кривої (Метод Ейлера)", marker='o', markersize=4, linestyle='-', color='blue')
    plt.title("Розв'язок y' = y*x^2 + 3*sin(y) (Метод Ейлера)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()
    plt.show()

elif choice == '2':
    k1_values, k2_values, k3_values, k4_values = [], [], [], []

    for i in range(1, len(x_values)):
        x_prev, y_prev = x_values[i - 1], y_values[-1]
        
        k1 = h * f(x_prev, y_prev)
        k2 = h * f(x_prev + 0.5 * h, y_prev + 0.5 * k1)
        k3 = h * f(x_prev + 0.5 * h, y_prev + 0.5 * k2)
        k4 = h * f(x_prev + h, y_prev + k3)
        
        # Метод Рунге-Кутта
        y_next = y_prev + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
        y_values.append(y_next)
        k1_values.append(k1)
        k2_values.append(k2)
        k3_values.append(k3)
        k4_values.append(k4)

    print(f"\nМетод Рунге-Кутта:\n{'i':<5}{'x_i':<10}{'y_i':<15}{'k1':<20}{'k2':<20}{'k3':<20}{'k4':<20}")
    print("-" * 130)
    for i in range(len(x_values)):
        x = x_values[i]
        y = y_values[i]
        k1 = k1_values[i - 1] if i > 0 else 0.0
        k2 = k2_values[i - 1] if i > 0 else 0.0
        k3 = k3_values[i - 1] if i > 0 else 0.0
        k4 = k4_values[i - 1] if i > 0 else 0.0
        
        print(f"{i:<5}{x:<10.2f}{y:<15.5f}{k1:<20.10f}{k2:<20.10f}{k3:<20.10f}{k4:<20.10f}")

    plt.figure(figsize=(10, 6))
    plt.plot(x_values, y_values, label="Наближений вигляд інтегральної кривої (Рунге-Кутта)", marker='o', markersize=4, linestyle='-', color='green')
    plt.title("Розв'язок y' = y*x^2 + 3*sin(y) (Метод Рунге-Кутта)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()
    plt.show()

else:
    print("Неправильний вибір. Введіть 1 або 2.")

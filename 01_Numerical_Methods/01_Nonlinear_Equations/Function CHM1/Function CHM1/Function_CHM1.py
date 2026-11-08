from sympy import symbols, diff, sin, log
import math
import sys

class Argument:
    def __init__(self, value):
        self.value = value
        self.n = self.get_n_input()
        self.first_significant_digit = self.first_significant_digit(self.value)
        self.rel_error = self.relative_error()
        self.abs_error = self.absolute_error()

    def first_significant_digit(self, x):
        if x == 0:
            raise ValueError("Не можна визначити першу значущу цифру для нуля.")
        
        x_abs = abs(x)
        if x_abs >= 1:
            while x_abs >= 10:
                x_abs /= 10
            return int(x_abs)
        else:
            x_str = str(x_abs)
            if '.' in x_str:
                x_str = x_str.lstrip('0').lstrip('.')

            x_str = x_str.lstrip('0')
            return int(x_str[0])

    def significant_digits_count(self, value):
        value_str = str(value)
        if '.' in value_str:
            before_decimal, after_decimal = value_str.split('.')
            before_decimal = before_decimal.lstrip('0')
            after_decimal = after_decimal.rstrip('0')
            return len(before_decimal) + len(after_decimal)
        else:
            return len(value_str.lstrip('0'))

    def get_n_input(self):
        print(f"Введіть кількість точних цифр для {self.value}: ")
        max_digits = self.significant_digits_count(self.value)
        print(f"!!!\nЗауважте: к-сть значущих цифр: {max_digits}")
        print(f"Тому 0<n<={max_digits}")
        while True:
            try:
                n = int(input())
                if n <= 0:
                    raise ValueError("Кількість точних цифр повинна бути більше 0.")
                if n > max_digits:
                    raise ValueError(f"Кількість точних цифр не може перевищувати {max_digits}.")
                return n
            except ValueError as e:
                print(e)

    def relative_error(self):
        return (1 / self.first_significant_digit) * (1 / 10)**(self.n - 1)

    def absolute_error(self):
        return self.value * self.rel_error


def compute_absolute_error(x1_val, x2_val, x3_val, abs_error_x1, abs_error_x2, abs_error_x3):

    dF_dx1_val = abs(dF_dx1.subs({x1: x1_val, x2: x2_val, x3: x3_val}))
    dF_dx2_val = abs(dF_dx2.subs({x1: x1_val, x2: x2_val, x3: x3_val}))
    dF_dx3_val = abs(dF_dx3.subs({x1: x1_val, x2: x2_val, x3: x3_val}))

    print(f"Часткові похідні функції в точці (x1={x1_val}, x2={x2_val}, x3={x3_val}):")
    print(f"∂F/∂x1 = {dF_dx1_val:.6f}")
    print(f"∂F/∂x2 = {dF_dx2_val:.6f}")
    print(f"∂F/∂x3 = {dF_dx3_val:.6f}")
    
    term1 = (dF_dx1_val * abs_error_x1)
    term2 = (dF_dx2_val * abs_error_x2)
    term3 = (dF_dx3_val * abs_error_x3)
    abs_error = term1 + term2 + term3

    print(f"Абсолютна похибка: {term1:.6f} (від x1) + {term2:.6f} (від x2) + {term3:.6f} (від x3) = {abs_error:.6f}")
    return abs_error

def compute_relative_error(x1_val, x2_val, x3_val, abs_error_x1, abs_error_x2, abs_error_x3):
   
    dlog_F_dx1_val = dlog_F_dx1.subs({x1: x1_val, x2: x2_val, x3: x3_val})
    dlog_F_dx2_val = dlog_F_dx2.subs({x1: x1_val, x2: x2_val, x3: x3_val})
    dlog_F_dx3_val = dlog_F_dx3.subs({x1: x1_val, x2: x2_val, x3: x3_val})

    print(f"Часткові похідні логарифма функції в точці (x1={x1_val}, x2={x2_val}, x3={x3_val}):")
    print(f"∂ln(F)/∂x1 = {dlog_F_dx1_val:.6f}")
    print(f"∂ln(F)/∂x2 = {dlog_F_dx2_val:.6f}")
    print(f"∂ln(F)/∂x3 = {dlog_F_dx3_val:.6f}")

    rel_error = (abs(dlog_F_dx1_val) * abs_error_x1 +
                 abs(dlog_F_dx2_val) * abs_error_x2 +
                 abs(dlog_F_dx3_val) * abs_error_x3)
    
    print(f"Відносна похибка: {rel_error * 100:.6f}%")
    return rel_error

print("=" * 40)
print("Група: ОІ-23")
print("Виконав: Гончаров Вадим Олегович")
print("Номер ЛР: №1")
print("Варіант: №4")
print("=" * 40)

x1, x2, x3 = symbols('x1 x2 x3')
F = 4*x1**2 + 5*x2**2 + 3*x3**2 - 4*x1*x2 - 2*x1 - sin(x1/x2)

dF_dx1 = diff(F, x1)
dF_dx2 = diff(F, x2)
dF_dx3 = diff(F, x3)

log_F = log(F)

dlog_F_dx1 = diff(log_F, x1)
dlog_F_dx2 = diff(log_F, x2)
dlog_F_dx3 = diff(log_F, x3)

print("Часткові похідні функції F:")
print(f"∂F/∂x1: {dF_dx1}")
print(f"∂F/∂x2: {dF_dx2}")
print(f"∂F/∂x3: {dF_dx3}")
print()

print("Часткові похідні логарифма функції F:")
print(f"∂ln(F)/∂x1: {dlog_F_dx1}")
print(f"∂ln(F)/∂x2: {dlog_F_dx2}")
print(f"∂ln(F)/∂x3: {dlog_F_dx3}")
print()

while True:
    option = input("**Оберіть умову виконання програми**\nA) Якщо користувач сам вводить к-сть точних цифр\nB) Якщо відома абсолютної похибка ( 𐤃=4*10^(-3) ): ").upper()
    
    if option == 'B':
        general_abs_error = 0.004

        x1_val = float(input("Введіть значення для x1: "))
        while True:
            x2_val = float(input("Введіть значення для x2 (не дорівнює нулю): "))
            if x2_val == 0:
                print("Помилка: x2 не може бути нульовим. Будь ласка, спробуйте ще раз.")
            else:
                break

        x3_val = float(input("Введіть значення для x3: "))

        abs_error_F = compute_absolute_error(x1_val, x2_val, x3_val, general_abs_error, general_abs_error, general_abs_error)
        print(f"Загальна абсолютна похибка функції: {abs_error_F:.10f}")
        rel_error_F = compute_relative_error(x1_val, x2_val, x3_val, general_abs_error, general_abs_error, general_abs_error)
        print(f"Відносна похибка функції: {rel_error_F*100:.6f}%")

        sys.exit()
    elif option == 'A':
        break
    else:
        print("Неправильний варіант. Будь ласка, введіть 'A' або 'B'.")

print()


while True:
    choice = input("**Оберіть кортеж аргументів**\n1) Ввести власні значення\n2) Використати значення за замовчуванням: ")
    
    if choice == '1':
        x1_val = float(input("Введіть значення для x1: "))
        while True:
            x2_val = float(input("Введіть значення для x2 (відмінне від нуля): "))
            if x2_val == 0:
                print("Помилка: x2 не може бути нульовим. Будь ласка, спробуйте ще раз.")
            else:
                break

        x3_val = float(input("Введіть значення для x3: "))
        break
    elif choice == '2':
        x1_val = 5.67321
        x2_val = 2.709457
        x3_val = 0.7811
        break
    else:
        print("Помилка. Будь ласка, введіть '1' або '2'.")

arg1 = Argument(x1_val)
arg2 = Argument(x2_val)
arg3 = Argument(x3_val)

print()
print("=" * 30)
print("Інформація щодо значення аргументів:")
print("=" * 30)

print(f"\nАргумент x1={arg1.value} | Відносна похибка: {arg1.rel_error * 100:.4f}% | Абсолютна похибка: {arg1.abs_error:.4f} | Перша значуща цифра: {arg1.first_significant_digit}")
print(f"Аргумент x2={arg2.value} | Відносна похибка: {arg2.rel_error * 100:.4f}% | Абсолютна похибка: {arg2.abs_error:.4f} | Перша значуща цифра: {arg2.first_significant_digit}")
print(f"Аргумент x3={arg3.value} | Відносна похибка: {arg3.rel_error * 100:.4f}% | Абсолютна похибка: {arg3.abs_error:.4f} | Перша значуща цифра: {arg3.first_significant_digit}")
print()
print("=" * 30)

abs_error_F = compute_absolute_error(x1_val, x2_val, x3_val, arg1.abs_error, arg2.abs_error, arg3.abs_error)
print(f"\nГранична абсолютна похибка функції: {abs_error_F:.10f}")

rel_error_F = compute_relative_error(x1_val, x2_val, x3_val, arg1.abs_error, arg2.abs_error, arg3.abs_error)
print(f"\nГранична відносна похибка функції F: {rel_error_F * 100:.6f}%")

F_val = F.subs({x1: x1_val, x2: x2_val, x3: x3_val})
print(f"\nГрачина відносна похибка функції (у виді 𐤃y/|F|): {(abs_error_F/F_val)*100:.6f}%")
print(f"\nЗначення функції F: {F_val}")
print("=" * 40)
print("Група: ОІ-23")
print("Виконав: Гончаров Вадим Олегович")
print("Номер ЛР: №1")
print("Варіант: №4")
print("=" * 40)
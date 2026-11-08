from sympy import symbols, diff, Poly
print("=" * 40)
print("Група: ОІ-23")
print("Виконав: Гончаров Вадим Олегович")
print("Номер ЛР: №2")
print("Варіант: №4")
print("=" * 40)
x = symbols('x')
polynomial = x**4 - 4*x**3 + 2*x - 8
derivative = diff(polynomial, x)
initial_degree = Poly(polynomial).degree()

def sturm_sequence(polynomial, derivative, degree):
    sequence = [Poly(polynomial), Poly(derivative)]
    while True:
        remainder = sequence[-2].div(sequence[-1])[1]
        degree -= 1
        if degree == 0:
            break
        sequence.append(remainder)
    return sequence

def sign_changes(sequence, value):
    signs = []
    for poly in sequence:
        signs.append(poly.eval(value))

    count_changes = 0

    for i in range(len(signs) - 1):
        if (signs[i] > 0) != (signs[i + 1] > 0):
            count_changes += 1

    return count_changes

def count_real_roots(polynomial, derivative, a, b, degree):
    sturm_seq = sturm_sequence(polynomial, derivative, initial_degree)
    return sign_changes(sturm_seq, a) - sign_changes(sturm_seq, b)

def upper_bound(polynomial):
    coeffs = Poly(polynomial).all_coeffs()
    an = abs(coeffs[0])

    max_coeff = 0 
    for i in range(len(coeffs)):
        if coeffs[i] < 0:
            max_coeff = max(max_coeff, abs(coeffs[i]))

    k = None
    for i in range(len(coeffs)):
        if coeffs[i] < 0:
            k = i
            break

    return (1 + (max_coeff / an) ** (1 / k))

def lower_bound(polynomial):
    neg_polynomial = polynomial.subs(x, -x)
    coeffs = Poly(neg_polynomial).all_coeffs()
    
    k = None
    for i in range(len(coeffs)):
        if coeffs[i] < 0:
            k = i
            break
    
    an = abs(coeffs[0])
    
    max_coeff = 0
    for i in range(len(coeffs)):
        if coeffs[i] < 0:
            max_coeff = max(max_coeff, abs(coeffs[i]))
    
    lower_bound_value = -(1 + (max_coeff / an) ** (1 / k))
    
    return lower_bound_value


upper_bound_value = upper_bound(polynomial)
print(f'Верхня границя для коренів: {upper_bound_value}')

lower_bound_value = lower_bound(polynomial)
print(f'Нижня границя для коренів: {lower_bound_value}')

num_real_roots = count_real_roots(polynomial, derivative, lower_bound_value, upper_bound_value, initial_degree)
print(f'Кількість дійсних коренів: {num_real_roots}')
print(f'Кількість комплексних коренів: {initial_degree-num_real_roots}')
print("=" * 40)
print("Група: ОІ-23")
print("Виконав: Гончаров Вадим Олегович")
print("Номер ЛР: №2")
print("Варіант: №4")
print("=" * 40)
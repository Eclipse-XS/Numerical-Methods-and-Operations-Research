def input_system():
    n = int(input("Enter the number of equations (greater than 2): "))
    while n <= 2:
        print("The number of equations must be greater than 2.")
        n = int(input("Enter the number of equations: "))
    
    A = []
    B = []

    print("Enter the coefficients of matrix A:")
    for i in range(n):
        row = []
        for j in range(n):
            row.append(float(input(f"A[{i+1},{j+1}] = ")))
        A.append(row)

    print("Enter the vector B:")
    for i in range(n):
        B.append(float(input(f"B[{i+1}] = ")))

    return A, B, n

def is_diagonally_dominant(A):
    n = len(A)
    for i in range(n):
        diagonal = abs(A[i][i])
        sum_off_diagonal = sum(abs(A[i][j]) for j in range(n) if i != j)
        if diagonal <= sum_off_diagonal:
            return False
    return True

def make_diagonally_dominant(A, B):
    n = len(A)
    for i in range(n):
        if abs(A[i][i]) <= sum(abs(A[i][j]) for j in range(n) if j != i):
            for k in range(i+1, n):
                if abs(A[k][i]) > sum(abs(A[k][j]) for j in range(n) if j != i):
                    A[i], A[k] = A[k], A[i]
                    B[i], B[k] = B[k], B[i]
                    break
    return A, B

def simple_iteration(A, B, epsilon, max_iterations=100):
    n = len(A)
    x = [0.0] * n  # Initial guess x = [0, 0, ..., 0]
    x_new = [0.0] * n

    for iteration in range(max_iterations):
        for i in range(n):
            sum_ax = sum(A[i][j] * x[j] for j in range(n) if i != j)
            x_new[i] = (B[i] - sum_ax) / A[i][i]

        norm = max(abs(x_new[i] - x[i]) for i in range(n))
        
        print(f"Iteration {iteration}: x = {x_new}, norm = {norm}")

        if norm < epsilon:
            print(f"Converged after {iteration + 1} iterations.")
            return x_new

        x = x_new[:]
    
    print(f"Did not converge after {max_iterations} iterations.")
    return x_new

def main():
    A, B, n = input_system()
    
    if not is_diagonally_dominant(A):
        print("Diagonal dominance condition is not satisfied, rearranging rows...")
        A, B = make_diagonally_dominant(A, B)
    
    epsilon = float(input("Enter the value of epsilon (less than 1): "))
    while epsilon >= 1:
        print("Epsilon must be less than 1.")
        epsilon = float(input("Enter the value of epsilon: "))

    result = simple_iteration(A, B, epsilon)
    print(f"Result: x = {result}")

if __name__ == "__main__":
    main()

import math

def is_diagonally_dominant(matrix):
    n = len(matrix)
    for i in range(n):
        diagonal_element = abs(matrix[i][i])
        sum_of_other_elements = sum(abs(matrix[i][j]) for j in range(n) if j != i)
        if diagonal_element <= sum_of_other_elements:
            return False
    return True

def vector_difference(first_vector, second_vector):
    n = len(first_vector)
    return_vector = [first_vector[i] - second_vector[i] for i in range(n)]
    return return_vector

def vector_first_norm(vector):
    return sum([abs(i) for i in vector])

def vector_second_norm(vector):
    return math.sqrt(sum([i * i for i in vector]))

def vector_max_norm(vector):
    return max([abs(i) for i in vector])

def correction_check(exact_solution, iter_solution):
    numer = vector_max_norm(vector_difference(exact_solution, iter_solution))
    denom = vector_max_norm(exact_solution)
    return numer / denom

def print_matrix(matrix, precision=0):
    for row in matrix:
        for col in row:
            spaces = 7
            template = '{:.' + str(precision) + 'f}'
            spaces -= (len(template.format(col)) - 1)
            print(template.format(col), end = " " * spaces)
        print()
    print()

def print_vector(vector, precision=0):
    for elem in vector:
        spaces = 7
        template = '{:.' + str(precision) + 'f}'
        spaces -= (len(template.format(elem)) - 1)
        print(template.format(elem), end = " " * spaces)
    print()

def print_report(solution_vector, iteration, method_name, is_max=False):
    print()
    print(method_name)
    print()
    if is_max:
        print("The answer did not reach the required tolerance when performing the maximum number of iterations.")
    print("Solution vector:")
    print_vector(solution_vector, 3)
    print("Iteration:\n", iteration, sep="")

def jacobi_method(A, b, initial_x, eps=1e-10, max_iterations=1000):
    '''
    Solving a system of linear equations using the Jacobi method.

    Matrix A is matrix of coefficients of the linear system.
    Vector b is right side of the linear system.
    Vector initial_x is initial guess of the linear system.
    Variable eps is tolerance.
    
    Function returns a vector that approximately solves a linear system.
    '''

    n = len(b)
    x = initial_x.copy()

    for k in range(max_iterations):
        x_old = x.copy()
        for i in range(n):
            sigma = 0
            for j in range(n):
                if j != i:
                    sigma += A[i][j] * x_old[j]
            x[i] = (b[i] - sigma) / A[i][i]
        if max(abs(x[i] - x_old[i]) for i in range(n)) < eps:
            print_report(x, k + 1, "Jacobi method")
            return x
    print_report(x, k + 1, "Jacobi method", True)
    return x

def gauss_seidel_method(A, b, initial_x, eps=1e-10, max_iterations=1000):
    '''
    Solving a system of linear equations using the Gauss-Seidel method.

    Matrix A is matrix of coefficients of the linear system.
    Vector b is right side of the linear system.
    Vector initial_x is initial guess of the linear system.
    Variable eps is tolerance.
    
    Function returns a vector that approximately solves a linear system.
    '''
    
    n = len(b)
    x = initial_x.copy()

    for k in range(max_iterations):
        x_old = x.copy()
        for i in range(n):
            left_sigma = 0
            right_sigma = 0
            for j in range(n):
                if j > i:
                    right_sigma += A[i][j] * x_old[j]
                elif j < i:
                    left_sigma += A[i][j] * x[j]
            x[i] = (b[i] - left_sigma - right_sigma) / A[i][i]
        if max(abs(x[i] - x_old[i]) for i in range(n)) < eps:
            print_report(x, k + 1, "Gauss-Seidel method")
            return x
    print_report(x, k + 1, "Gauss-Seidel method", True)
    return x

def successive_over_relaxation(A, b, initial_x, w, eps=1e-10, max_iterations=1000):
    '''
    Solving a system of linear equations using the successive over-relaxation.

    Matrix A is matrix of coefficients of the linear system.
    Vector b is right side of the linear system.
    Vector initial_x is initial guess of the linear system.
    
    Function returns a vector that approximately solves a linear system.
    '''
    
    n = len(b)
    x = initial_x.copy()

    for k in range(max_iterations):
        x_old = x.copy()
        for i in range(n):
            left_sigma = 0
            right_sigma = 0
            for j in range(n):
                if j > i:
                    right_sigma += A[i][j] * x_old[j]
                elif j < i:
                    left_sigma += A[i][j] * x[j]
            x[i] = (1 - w) * x_old[i] + (w * (b[i] - left_sigma - right_sigma) / A[i][i])
        if max(abs(x[i] - x_old[i]) for i in range(n)) < eps:
            print_report(x, k + 1, "Successive over-relaxation with w = " + str(w))
            return x
    print_report(x, k + 1, "Successive over-relaxation with w = " + str(w), True)
    return x

n = 15

A = [[0 for i in range(n)] for i in range(n)]
b = [0 for i in range(n)]
x_initial = b.copy()
x_exact = [i + 1 for i in range(n)]

for i in range(n):
    for j in range(n):
        A[i][j] = 0.01 * (math.sqrt(i + 1) + math.sqrt(j + 1))
for i in range(n):
    A[i][i] = 5 * math.sqrt(i + 1)

for i in range(n):
    for j in range(n):
        b[i] += (A[i][j] * x_exact[j])

precision = 3
max_iterations = 1000
tolerance = 1e-10

print_matrix(A, precision)
print_vector(x_exact, precision)
print_vector(b, precision)

x_jacobi = jacobi_method(A, b, x_initial, tolerance, max_iterations)
print("Check report:", correction_check(x_exact, x_jacobi))
x_gauss_seidel = gauss_seidel_method(A, b, x_initial, tolerance, max_iterations)
print("Check report:", correction_check(x_exact, x_gauss_seidel))
x_successive_over_relaxation_0_5 = successive_over_relaxation(A, b, x_initial, 0.5, tolerance, max_iterations)
print("Check report:", correction_check(x_exact, x_successive_over_relaxation_0_5))
x_successive_over_relaxation_1_5 = successive_over_relaxation(A, b, x_initial, 1.5, tolerance, max_iterations)
print("Check report:", correction_check(x_exact, x_successive_over_relaxation_1_5))
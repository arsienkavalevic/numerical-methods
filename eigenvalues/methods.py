import numpy as np

def generate_matrix(dim):
    random_matrix = np.random.random_integers(-10,10,size=(dim,dim))
    a = (random_matrix + random_matrix.T) // 2
    return a

def jacobi_rotation_method(A, eps, max_iterations=1000):
    n = len(A)
    iter_num = 0

    P = np.eye(n)

    for k in range(max_iterations):
        max_off_diag = 0.0
        max_i, max_j = 0, 0
        for i in range(n):
            for j in range(i + 1, n):
                if np.abs(A[i, j]) > max_off_diag:
                    max_off_diag = np.abs(A[i, j])
                    max_i, max_j = i, j
        if max_off_diag < eps:
            break

        theta = 0.5 * np.arctan2(2 * A[max_i, max_j], A[max_i, max_i] - A[max_j, max_j])

        c = np.cos(theta)
        s = np.sin(theta)
        V = np.eye(n)
        V[max_i, max_i] = c
        V[max_j, max_j] = c
        V[max_i, max_j] = -s
        V[max_j, max_i] = s

        A = np.dot(np.dot(V.T, A), V)
        P = np.dot(P, V)

        iter_num = k + 1

    eigenvalues = np.diag(A)
    eigenvectors = P.T

    return iter_num, eigenvectors, eigenvalues

def power_method(A, eps, max_iterations=1000):
    n = len(A)
    iter_num = 0
    eigenvector = np.random.randint(-10, 10, size=n)
    eigenvector = eigenvector - np.linalg.norm(eigenvector)
    for k in range(max_iterations):
        y = np.dot(A, eigenvector)
        eigenvector_new = y / np.linalg.norm(y)
        eigenvalue = np.dot(eigenvector, np.dot(A, eigenvector)) / np.dot(eigenvector, eigenvector)

        if np.linalg.norm(np.dot(A, eigenvector) - np.dot(eigenvalue, eigenvector)) < eps:
            break

        eigenvector = eigenvector_new
        iter_num = k + 1

    return iter_num, eigenvector, eigenvalue


DIM = 3
EPS = 1e-6
A = generate_matrix(DIM)
print("Matrix: ")
print(A)
jacobi_iter_num, jacobi_eigenvectors, jacobi_eigenvalues = jacobi_rotation_method(A, EPS)
print("\nNum of iterations with jacobi rotation method: ", end="")
print(jacobi_iter_num)
print("\nEigenvectors with jacobi rotation method: ")
print(jacobi_eigenvectors)
print("\nEigenvalues with jacobi rotation method: ")
print(jacobi_eigenvalues)
print("\nJacobi method norm errors: ")
for i in range(DIM):
    print(np.linalg.norm((np.dot(A, jacobi_eigenvectors[i]) - np.dot(jacobi_eigenvalues[i], jacobi_eigenvectors[i]))))

power_iter_num, power_eigenvector, power_eigenvalue = power_method(A, EPS)
print("\nNum of iterations with power method: ", end="")
print(power_iter_num)
print("\nEigenvector with power method: ")
print(power_eigenvector)
print("\nEigenvalue with power method: ")
print(power_eigenvalue)
print("\nPower method norm error: ")
print(np.linalg.norm((np.dot(A, power_eigenvector) - np.dot(power_eigenvalue, power_eigenvector))))
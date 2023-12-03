import math

def summation(left, right, i, j, bottom, top):
    rtn = 0
    for k in range(bottom, top + 1):
        rtn += (left[i][k] * right[k][j])
    return rtn

def determinant(l, u):
    n = len(l)
    l_det = 1
    for i in range(n):
        l_det *= l[i][i]
    u_det = 1
    for i in range(n):
        u_det *= u[i][i]
    return l_det * u_det

def matrix_invertation(l, u):
    n = len(l)
    d = [[0 for i in range(n)] for j in range(n)]
    e = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        e[i][i] = 1
    for i in range(n):
        d[0][i] = e[0][i]
    for j in range(n):
        for i in range(1, n):
            d[i][j] = e[i][j] - summation(l, d, i, j, 0, i - 1)
    x = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        x[n - 1][i] = d[n - 1][i] / u[n - 1][n - 1]
    for j in range(n - 1, -1, -1):
        for i in range(n - 2, -1, -1):
            x[i][j] = (d[i][j] - summation(u, x, i, j, i + 1, n - 1)) / u[i][i]
    return x

def matrix_multiplication(matrix1, matrix2):
    n = len(matrix1)
    result_matrix = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            result_matrix[i][j] = sum([matrix1[i][k] * matrix2[k][j] for k in range(n)])
    return result_matrix

def matrix_summation(matrix1, matrix2):
    n = len(matrix1)
    result_matrix = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            result_matrix[i][j] = matrix1[i][j] + matrix2[i][j]
    return result_matrix

def matrix_subtraction(matrix1, matrix2):
    n = len(matrix1)
    result_matrix = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            result_matrix[i][j] = matrix1[i][j] - matrix2[i][j]
    return result_matrix

def matrix_residual_norm(matrix1, matrix2):
    matrix_residual = matrix_subtraction(matrix1, matrix2)
    columns = []
    n = len(matrix_residual)
    for j in range(n):
        sum = 0
        for i in range(n):
            sum += matrix[i][j]
        columns.append(sum)
    return max(columns)

matrix = [[10,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0],
          [ 1, 10,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0],
          [ 0,  1, 10,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0],
          [ 0,  0,  1, 10,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0],
          [ 0,  0,  0,  1, 10,  1,  0,  0,  0,  0,  0,  0,  0,  0, 0],
          [ 0,  0,  0,  0,  1, 10,  1,  0,  0,  0,  0,  0,  0,  0, 0],
          [ 0,  0,  0,  0,  0,  1, 10,  1,  0,  0,  0,  0,  0,  0, 0],
          [ 0,  0,  0,  0,  0,  0,  1, 10,  1,  0,  0,  0,  0,  0, 0],
          [ 0,  0,  0,  0,  0,  0,  0,  1, 10,  1,  0,  0,  0,  0, 0],
          [ 0,  0,  0,  0,  0,  0,  0,  0,  1, 10,  1,  0,  0,  0, 0],
          [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1, 10,  1,  0,  0, 0],
          [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, 10,  1,  0, 0],
          [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, 10,  1, 0],
          [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, 10, 1],
          [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 1]]
matrix_len = len(matrix)

d_array = list(range(1, 16))

#A matrix output
print("A =")
for row in matrix:
    for col in row:
        spaces = 2
        template = '{:.' + str(0) + 'f}'
        spaces -= (len(template.format(col)) - 1)
        print(template.format(col), end = " " * spaces)
    print()
print()

#b vector output
print("b =")
print(d_array)
print()

#-------------------------------------------#
#-first-part-of-ridiagonal-matrix-algorithm-#
#-------------------------------------------#

a_array = [0]
for i in range(1, len(matrix[0])):
    a_array.append(matrix[i][i - 1])

b_array = []
for i in range(len(matrix[0])):
    b_array.append(matrix[i][i])

c_array = []
for i in range(1, len(matrix[0])):
    c_array.append(matrix[i - 1][i])
c_array.append(0)

xi_array = []
eta_array = []

xi_array.append(c_array[0] / (-1 * b_array[0]))

for i in range(1, len(matrix[0])):
    xi = (c_array[i] / ((-1 * b_array[i]) - (a_array[i] * xi_array[i - 1])))
    xi_array.append(xi)

eta_array.append(d_array[0] / b_array[0])

for i in range(1, len(matrix[0])):
    eta = (((a_array[i] * eta_array[i - 1]) - d_array[i]) / ((-1 * b_array[i]) - (a_array[i] * xi_array[i - 1])))
    eta_array.append(eta)

#--------------------------------------------#
#-second-part-of-ridiagonal-matrix-algorithm-#
#--------------------------------------------#

xi_array = xi_array[::-1]
eta_array = eta_array[::-1]

y_array = []

y_array.append(eta_array[0])

for i in range(1, len(matrix[0])):
    y = (xi_array[i] * y_array[i - 1]) + eta_array[i]
    y_array.append(y)

#y vector output
print("y =")
print(y_array)
print()

#------------------#
#-LU-decomposition-#
#------------------#

l_matrix = [[0 for i in range(len(matrix))] for j in range(len(matrix))]
u_matrix = [[0 for i in range(len(matrix))] for j in range(len(matrix))]

for i in range(len(d_array)):
    for j in range(len(d_array)):
        u_matrix[i][j] = 0
        l_matrix[i][j] = 0
    l_matrix[i][i] = 1

for i in range(len(d_array)):
    for j in range(len(d_array)):
        if i <= j:
            u_matrix[i][j] = matrix[i][j] - summation(l_matrix, u_matrix, i, j, 0, i)
        if i > j:
            l_matrix[i][j] = (matrix[i][j] - summation(l_matrix, u_matrix, i, j, 0, j)) / u_matrix[j][j]

#L matrix output
print("L =")
for row in l_matrix:
    for col in row:
        spaces = 7
        template = '{:.' + str(4) + 'f}'
        spaces -= (len(template.format(col)) - 1)
        print(template.format(col), end = " " * spaces)
    print()
print()

#U matrix output
print("U =")
for row in u_matrix:
    for col in row:
        spaces = 7
        template = '{:.' + str(4) + 'f}'
        spaces -= (len(template.format(col)) - 1)
        print(template.format(col), end = " " * spaces)
    print()
print()

#matrix_invertation method call
invert_matrix = matrix_invertation(l_matrix, u_matrix)

#inverted A matrix output
print("A^(-1) =")
for row in invert_matrix:
    for col in row:
        spaces = 7 
        template = '{:.' + str(4) + 'f}'
        spaces -= (len(template.format(col)) - 1)
        print(template.format(col), end = " " * spaces)
    print()
print()

#multiplication of A matrix and inverted A matrix output
print("A * A^(-1) = E =")
for row in matrix_multiplication(matrix, matrix_invertation(l_matrix, u_matrix)):
    for col in row:
        spaces = 2
        template = '{:.' + str(0) + 'f}'
        spaces -= (len(str(template.format(abs(col)))) - 1)
        print(template.format(abs(col)), end = " " * spaces)
    print()
print()

zero_matrix = [[0 for i in range(matrix_len)] for j in range(matrix_len)]
e_matrix = [[0 for i in range(matrix_len)] for j in range(matrix_len)] #identity matrix
for i in range(matrix_len):
    e_matrix[i][i] = 1

#LU - A norm output
print("LU - A norm:")
print(matrix_residual_norm(matrix_subtraction(matrix_multiplication(l_matrix, u_matrix), matrix), zero_matrix))
print()

#Residual norm output
print("Residual norm:")
print(matrix_residual_norm(matrix_multiplication(matrix, invert_matrix), e_matrix))
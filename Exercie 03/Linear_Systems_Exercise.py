import numpy as np
import matplotlib.pyplot as plt
m = [[1,2,3],[4,5,6],[7,8,9]]
n = [[10],[11],[12]]
def mat_mul(m, n):
    if len(m[0]) != len(n): return print("Dimensions do not match for Matrix Multiplication")
    Matrix = [[0 for _ in range(len(n[0]))] for _ in range(len(m))]
    Matrix = [[sum(m[i][k] * n[k][j] for k in range(len(n)) for j in range(len(n[0])))] for i in range(len(m))]
    return Matrix
Matrix = mat_mul(m,n)
print(Matrix)


def get_index(i, j, N):
    return i * N + j

def back_sub(M, b, N):
    x = np.zeros(N)
    for i in range(N - 1, -1, -1):
        x[i] = b[i]
        for j in range(i + 1, N):
            x[i] -= M[get_index(i, j, N)] * x[j]
        
        x[i] /= M[get_index(i, i, N)]
    
    return x

N = 3
M = np.array([1, 2, 3, 0, 5, 6, 0, 0, 9], dtype=float)
b = np.array([10, 11, 12], dtype=float)
x = back_sub(M, b, N)
print("Solution vector x:", x)
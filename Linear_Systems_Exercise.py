import numpy as np
import matplotlib.pyplot as plt
m = [[1,2],[5,6],[8,9]]
n = [[10],[11],[12]]
def mat_mul(m, n):
    if len(m[0]) != len(n): return print("Dimensions do not match for Matrix Multiplication")
    Matrix = [[0 for _ in range(len(n[0]))] for _ in range(len(m))]
    Matrix = [[sum(m[i][k] * n[k][j] for k in range(len(n)) for j in range(len(n[0])))] for i in range(len(m))]
    return Matrix
Matrix = mat_mul(m,n)
print(Matrix)
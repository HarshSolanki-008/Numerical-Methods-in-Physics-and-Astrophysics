import numpy as np

def random_matrix(m,n):
    A = np.zeros((m,n))
    B = np.zeros(n)
    for i in range(m):
        for j in range(n):
            A[i][j] = np.random.rand()*10
            B[j] = np.random.rand()*10
    return A,B

m,n = 3,3
A,B = random_matrix(m,n)

print(A)
print(B)
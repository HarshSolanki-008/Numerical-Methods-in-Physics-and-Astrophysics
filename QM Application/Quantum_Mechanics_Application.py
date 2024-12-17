import numpy as np

def array_fill(fn,a,b,n):
    x = np.linspace(a,b,n)
    y = fn(x)
    return x,y

def Newton_Cootes_2nd_Order(x,y,a,b,n):
    m = n // 2
    if n % 2 == 0:
        print("Must be odd number of points!")
        exit
    else:
        h = (x[-1] - x[0]) / (n - 1)
        solution = y[0] + y[-1] + 4 * np.sum(y[1:n-1:2]) + 2 * np.sum(y[2:n-2:2])
        solution = solution * (h / 3)
        fourth_derivative =fn(x[m + 2]) - 4*fn(x[m + 1]) + 6*fn(x[m]) - 4*fn(x[m - 1]) + fn(x[m - 2])
        error = (1/90)*fourth_derivative
        return solution,error
    
def fn(x):
    return x + x**3
    
a = 0
b = 1
n = 5
x,y = array_fill(fn,a,b,n)
solution, error = Newton_Cootes_2nd_Order(x,y,a,b,n)
print("Solution is: ", solution)
print("Error is: ", error)

# def psi(x,t,L,omega1,omega2):
#    if 0 < x < L:
#        psi = (1/np.sqrt(L))*(np.sin(np.pi*x/L)*np.exp(-1j*omega1*t) + np.sin(2*np.pi*x/L)*np.exp(-1j*omega2*t))
#    else:
#        psi = 0
#    P = np.abs(psi)**2  
#    return P

def psi(x,t,L,omega1,omega2):
    psi = (1/np.sqrt(L))*(np.sin(np.pi*x/L)*np.exp(-1j*omega1*t) + np.sin(2*np.pi*x/L)*np.exp(-1j*omega2*t))
    return np.abs(psi)**2

L = 2
omega1 = 3
omega2 = 4.5
del_omega = omega2 - omega1

a = 3*L/4
b = L
n = 255
time = [0,np.pi/del_omega]

for t in time:
    x,y = array_fill(lambda x: psi(x,t,L,omega1,omega2),a,b,n)
    solution1,error1 = Newton_Cootes_2nd_Order(x,y,a,b,n)
    print(f"Solution for t = {t} is: {solution1}")
    print(f"Error here is: {error1}")
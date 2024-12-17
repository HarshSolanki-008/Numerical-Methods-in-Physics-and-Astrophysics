import numpy as np
import matplotlib.pyplot as plt

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


def psi(x,t,L,omega1,omega2):
    psi = (1/np.sqrt(L))*(np.sin(np.pi*x/L)*np.exp(-1j*omega1*t) + np.sin(2*np.pi*x/L)*np.exp(-1j*omega2*t))
    return np.abs(psi)**2

L = 2
omega1 = 3
omega2 = 4.5
del_omega = omega2 - omega1

a = 3*L/4
b = L
n = np.arange(5,503,2)
time = [0,np.pi/del_omega]

solutions_t1 = []
errors_t1 = []
h_t1 = []

solutions_t2 = []
errors_t2 = []
h_t2 = []

for i in range(len(n)):
    x,y = array_fill(lambda x: psi(x,time[0],L,omega1,omega2),a,b,n[i])
    solution1,error1 = Newton_Cootes_2nd_Order(x,y,a,b,n[i])
    solutions_t1.append(solution1)
    errors_t1.append(error1)
    h_t1.append(L/(n[i] - 1))

for i in range(len(n)):
    x,y = array_fill(lambda x: psi(x,time[1],L,omega1,omega2),a,b,n[i])
    solution1,error1 = Newton_Cootes_2nd_Order(x,y,a,b,n[i])
    solutions_t2.append(solution1)
    errors_t2.append(error1)
    h_t2.append(L/(n[i] - 1))

plt.figure("log(E) vs log(h)")
plt.subplot(1,2,1)
plt.plot(np.log(errors_t1),np.log(h_t1), label = f"for time = {time[0]}", color = "red")
plt.xlabel("log(h)")
plt.ylabel("log(E)")
plt.legend()

plt.subplot(1,2,2)
plt.plot(np.log(errors_t2),np.log(h_t2), label = f"for time = {time[1]}", color = "black")
plt.xlabel("log(h)")
plt.ylabel("log(E)")
plt.legend()

plt.show()
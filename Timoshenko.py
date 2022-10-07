import numpy as np
import matplotlib.pyplot as plt

sigma_yield = 400e+9 #Pa
E = 200e+9 #Pa
nu = 0.3
rho = 7800 #kg/m^3
g = 9.81 #m/s^2
m = 100 #kg
h = 0.05 #m
b = 0.05 #m
G = E / (2*(1+nu)) #Pa
K = 5/6

L = [0.25, 2] #m

P = m * g #N
A = b * h #m^2
I = b * h**3/12 #m^4

def w(x, i):
    first = rho * A * g * (1 / (G * K * A) * (x**2 / 2 - L[i] * x) + 1 / (E * I) * (-x**4 / 24 + L[i] * x**3 / 6 - L[i]**2 * x**2 / 4))
    second = -P * (1 / (G * K * A) * x + 1 / (E * I) * (-x**3 / 6 + L[i] * x**2/2))
    return first + second
    ##########

def wEB(x, i):
    first = - rho * A * g / (E * I) * (x**4 / 24 - L[i] * x**3 / 6 + L[i]**2 * x**2/4)
    second = P / (E * I) * (x**3 / 6 - L[i] * x**2/ 2)
    return first + second

fig1,ax1 = plt.subplots()
x = np.linspace(0, L[0], 1000)
plt.plot(x, w(x, 0) * 1000, label = 'Timoshenko model') #mm
plt.plot(x, wEB(x, 0) * 1000, '--r', label = 'Euler-Bernoulli model') #mm
plt.title('Timoshenko and EB  beam deflection, L = 0.25m')
plt.legend(loc = 3)
plt.ylabel('w [mm]')
plt.xlabel('x [m]')

fig2,ax2 = plt.subplots()
plt.plot(x, w(x, 1) * 1000, label = 'Timoshenko model') #mm
plt.plot(x, wEB(x, 1) * 1000, '--r', label = 'Euler-Bernoulli model') #mm
plt.title('Timoshenko and EB beam deflection, L = 2m')
plt.legend(loc = 3)
plt.ylabel('w [mm]')
plt.xlabel('x [m]')

print(w(L[0], 0) * 1000)
print(w(L[1], 1) * 1000)
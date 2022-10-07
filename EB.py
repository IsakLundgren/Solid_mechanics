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

L = [0.25, 2] #m

P = m * g #N
A = b * h #m^2
I = b * h**3/12

def w(x, i):
    first = - rho * A * g / (E * I) * (x**4 / 24 - L[i] * x**3 / 6 + L[i]**2 * x**2/4)
    second = P / (E * I) * (x**3 / 6 - L[i] * x**2/ 2)
    return first + second

def sigma_xx(x, i): #at z = -h/2
    first = - rho * A * g * h / (2 * I) * (x**2 / 2 - L[i] * x + L[i]**2 / 2)
    second = - P * h / (2 * I) * (L[i] - x)
    return first + second

def sigma_eM(sigma_xx):
    return np.abs(sigma_xx^2)

fig1,ax1 = plt.subplots()
x = np.linspace(0, L[0], 1000)
plt.plot(x, w(x, 0) * 1000) #mm
plt.title('Euler Bernoulli beam deflection, L = 0.25m')
plt.ylabel('w [mm]')
plt.xlabel('x [m]')

fig2,ax2 = plt.subplots()
x = np.linspace(0, L[1], 1000)
plt.plot(x, w(x, 1) * 1000) #mm
plt.title('Euler Bernoulli beam deflection, L = 2m')
plt.ylabel('w [mm]')
plt.xlabel('x [m]')

fig3,ax3 = plt.subplots()
x = np.linspace(0, L[0], 1000)
plt.plot(x, sigma_xx(x, 0)) #Pa
plt.title('Normal stress at $z = -h/2$, L = 0.25m')
plt.ylabel('$\sigma_{xx}$ [mm]')
plt.xlabel('x [m]')

fig4,ax4 = plt.subplots()
x = np.linspace(0, L[1], 1000)
plt.plot(x, sigma_xx(x, 1)) #Pa
plt.title('Normal stress at $z = -h/2$, L = 2m')
plt.ylabel('$\sigma_{xx}$ [mm]')
plt.xlabel('x [m]')

print(w(L[0], 0) * 1000)
print(w(L[1], 1) * 1000)
print(sigma_xx(0, 0))
print(sigma_xx(0, 1))
plt.show()

#Gaming
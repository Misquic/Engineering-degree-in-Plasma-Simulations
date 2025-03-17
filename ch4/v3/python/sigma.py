import numpy as np
import numba as nb
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

s = 10
# Funkcja modelu: y = a * ln(x / b) / x
def model(x, a, b, c):
    return a * np.log(x / b) / x*np.exp(-c/x)

# Wczytanie danych z pliku
filename = "data/Oxygen_ionization.txt"  # Zamień na nazwę swojego pliku

data = np.loadtxt(filename)

# Rozdzielenie na x i y
x_data, y_data = data[:, 0], data[:, 1]

# Dopasowanie krzywej
initial_guess = [1e-20, 1e2, 1.1]  # Początkowe wartości dla a i b
params, covariance = curve_fit(model, x_data, y_data, p0=initial_guess)
errors = np.sqrt(np.diag(covariance))

# Wyodrębnienie parametrów
a, b, c = params
print(f"Dopasowane parametry: a = {a:.3e}, b = {b:.3e}, c = {c:.3e}")
print(errors)

# Obliczenie dopasowanej funkcji
x_fit = np.logspace(np.log10(min(x_data)/10), np.log10(max(x_data)*1000), 800)

y_fit = model(x_fit, a, b, c)

# Wizualizacja wyników
fig, ax = plt.subplots()
fig.set_size_inches(6,3)

ax.scatter(x_data, y_data, color='blue', label=r'$\sigma_{jon}$', s = s)
ax.plot(x_fit, y_fit, color='red', label=r"Dopasowana krzywa")#: $y = "+f"{a*1e18:.3}" + r"^{-18}\ln(\frac{x}{" + f"{b:.3}" + r"}) / x \cdot\exp(-\frac{" + f"{c:.3}" + r"}{x})$")
ax.set_xlabel('$E$ [eV]')
ax.set_ylabel(r'$\sigma$ [$\text{m}^2$]')
ax.grid()

# # Wielomian Czebyszewa
# def chebyshev_poly(n, x):

#     result = np.ones_like(x)
#     for i in range(len(x)):
#         el = x[i]
#         if el > 1:
#             result[i] = ((el + np.sqrt(el*el - 1))**n + (el - np.sqrt(el*el - 1))**n)*0.5

#         else:
#             result[i] = np.cos(n*np.arccos(el))
#     return result

# # Charakterystyka amplitudowa filtru Czebyszewa
# def chebyshev_filter(f, f_c, epsilon, n, a):
#     print(f"f_c: {f_c}, epsilon: {epsilon}, n: {n}, a: {a}")
#     x = f / f_c
#     Tn = chebyshev_poly(n, x)
#     return a / np.sqrt(1 + epsilon**2 * Tn**2)

# # Wczytanie danych z pliku
filename = "data/Oxygen_momentum_transfer.txt"  # Zamień na nazwę swojego pliku

data = np.loadtxt(filename)

# Rozdzielenie na x i y
x_data, y_data = data[:, 0], data[:, 1]

# initial_guess = [1.2e1, 5, 1.5, 5e-20]
# bounds = ([1e0, 1.0, 1.0, 1e-21],  # Dolne granice
#           [5e2, 50, 2.5, 6e-20]) 
# # Dopasowanie funkcji
# popt, pcov = curve_fit(lambda f, f_c, epsilon, n, a: chebyshev_filter(f, f_c, epsilon, n, a), 
#                        x_data, y_data, p0=initial_guess, bounds = bounds, method='trf')

# # Wyniki
# f_c, epsilon, n, a = popt
# print(f"Dopasowane parametry: f_c = {f_c}, epsilon = {epsilon}, n = {n}, a = {a}")

# x_fit = np.logspace(np.log10(min(x_data)), np.log10(max(x_data)), 500)
# # y_fit = chebyshev_filter(x_fit, *initial_guess)
# y_fit = chebyshev_filter(x_fit, *popt)

# # Wizualizacja wyników
# plt.figure(figsize=(8, 6))
ax.scatter(x_data, y_data, color='green', label=r'$\sigma_{całk}$', s = s)
# plt.plot(x_fit, y_fit, color='red', label=f'Dopasowana krzywa')
# plt.xlabel('x')
# plt.ylabel('y')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(y_data.min()*0.5, y_data.max()*1.5)
# plt.title('Dopasowanie krzywej')
# plt.legend()
# plt.grid()
# plt.show()
ax.legend()
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)

fig.savefig("data/sigmas.pdf", dpi = 350)
fig.savefig("data/sigmas.png", dpi = 350)

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Funkcja modelu: y = a * ln(x / b) / x
def model(x, a, b, c):
    return a * np.log(x / b) / x*np.exp(-c/x)

# Wczytanie danych z pliku
filename = "sigma.dat"  # Zamień na nazwę swojego pliku
data = np.loadtxt(filename)

# Rozdzielenie na x i y
x_data, y_data = data[:, 0], data[:, 1]

# Dopasowanie krzywej
initial_guess = [1e-20, 1e2, 1.1]  # Początkowe wartości dla a i b
params, covariance = curve_fit(model, x_data, y_data, p0=initial_guess)

# Wyodrębnienie parametrów
a, b, c = params
print(f"Dopasowane parametry: a = {a:.3e}, b = {b:.3e}, c = {c:.3e}")

# Obliczenie dopasowanej funkcji
x_fit = np.linspace(min(x_data), max(x_data), 500)
y_fit = model(x_fit, a, b, c)

# Wizualizacja wyników
plt.figure(figsize=(8, 6))
plt.scatter(x_data, y_data, color='blue', label='Dane eksperymentalne')
plt.plot(x_fit, y_fit, color='red', label=f'Dopasowana krzywa: $y = {a:.3e} \ln(x/{b:.3e}) / x$')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Dopasowanie krzywej')
plt.legend()
plt.grid()
plt.show()
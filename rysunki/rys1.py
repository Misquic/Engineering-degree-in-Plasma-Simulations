import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# Generowanie danych
# np.random.seed(45)
data = np.random.normal(loc=0, scale=400, size=100)  # Średnia=0, odchylenie standardowe=400
x = np.linspace(-1200, 1200, 500)

fig, ax = plt.subplots()
fig.set_size_inches(15,8)

# Histogram
ax.hist(data, bins=12, density=True, alpha=0.5, color='gray', edgecolor='black')

# Krzywa normalnego rozkładu
mean, std = np.mean(data), np.std(data)
pdf = norm.pdf(x, mean, std)
ax.plot(x, pdf, 'k-', linewidth=2)

# Punkty na osi
ax.scatter(data, np.zeros_like(data)+0.05*pdf.max(), facecolors='none', edgecolors='black', s = 50)

# Opis osi
ax.set_xlabel(r'prędkość $\frac{\text{m}}{\text{s}}$')
# ax.spines['left'].set_visible(False)
ax.yaxis.set_ticklabels([])
ax.yaxis.set_ticks([])

# Wyświetlanie wykresu
fig.savefig("obrazki2/fluid_vs_kinetic.png", dpi = 350)
fig.savefig("obrazki2/fluid_vs_kinetic.pdf", dpi = 350)
plt.show()


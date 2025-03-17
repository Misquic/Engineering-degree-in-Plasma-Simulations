import matplotlib.pyplot as plt
import re

# Wczytaj dane z pliku
def parse_data(filename):
    data_map = []
    data_vector = []
    
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith("used map") or line.startswith("used vector"):
                match = re.search(r"(used map|used vector).*frac: ([\d\.e-]+).*size ele: (\d+).*size neus: (\d+).*time: ([\d\.e-]+)", line)
                if match:
                    method, frac, size_ele, size_neus, time = match.groups()
                    frac = float(frac)
                    time = float(time)
                    
                    if "map" in method:
                        data_map.append((frac, time))
                    else:
                        data_vector.append((frac, time))
    
    return data_map, data_vector

# Rysowanie wykresu
def plot_data(data_map, data_vector):
    plt.figure(figsize=(10, 5))
    
    if data_map:
        fracs_map, times_map = zip(*data_map)
        plt.plot(fracs_map, times_map, 'b.', label='Used pointers')
    
    if data_vector:
        fracs_vector, times_vector = zip(*data_vector)
        plt.plot(fracs_vector, times_vector, 'r.', label='Used Vector')
    
    plt.xlabel('Frac')
    plt.ylabel('Time (s)')
    plt.title('Time vs Frac for Map and Vector')
    plt.legend()
    plt.grid()
    plt.show()

# Przykładowe użycie
filename = "outputs/output.txt"  # Podaj nazwę pliku z danymi
map_data, vector_data = parse_data(filename)
plot_data(map_data, vector_data)

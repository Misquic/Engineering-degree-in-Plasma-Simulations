# import numpy as np
# import matplotlib.pyplot as plt
# import vtk
# from vtk.util.numpy_support import vtk_to_numpy
# import sys

# if(len(sys.argv) > 1): 
#     scale = int(sys.argv[1])
# else:
#     scale = 1

# plt.rcParams["font.size"] = 14

# # Wczytaj plik .vtp
# reader = vtk.vtkXMLImageDataReader()
# num = str(f"{0:05}")
# # folder = "results_langmuir_2"
# # num_files = int(1500/scale)

# folder = "results"
# num_files = int(400/scale)




# reader.SetFileName(r"../" + folder +r"/fields_"+num+".vti")
# reader.Update()
# mesh = reader.GetOutput()
# dimensions = mesh.GetDimensions()
# point_data = mesh.GetPointData()

# for i in range(point_data.GetNumberOfArrays()):
#     if(point_data.GetArrayName(i) == "nd.e-"):
#         nd_e_index = i

# # num_files = 3004+1
# file_step = 30*scale
# nd_e_all = np.zeros((num_files, dimensions[0]))
# time = np.zeros(num_files)
# dt = 20e-9
# print(nd_e_all.shape)
# print(time.shape)

# if(len(sys.argv)>1):
#     for files_index in range(num_files):

#         reader = vtk.vtkXMLImageDataReader()
#         num = str(f"{files_index*file_step:05}")
#         print(num)

#         reader.SetFileName(r"../" + folder +r"/fields_"+num+".vti")
#         reader.Update()

#         # Pobierz dane
#         mesh = reader.GetOutput()
#         # print(mesh)
#         # print(f"mesh: {mesh}")
#         # print(f"dimensions: {dimensions}")

#         point_data = mesh.GetPointData()
#         nd_e = point_data.GetArray(nd_e_index)
#         # print(f"nd.e-: {nd_e}")
#         nd_e_numpy = vtk_to_numpy(nd_e).reshape(dimensions)
#         # print(nd_e_numpy)

#         nd_e_x = np.zeros(dimensions[0])
#         for k in range(dimensions[2]):
#             nd_e_x[k] += nd_e_numpy[int(dimensions[0]/2)][int(dimensions[1]/2)][k]
#         nd_e_x /= (dimensions[1]*dimensions[2])
#         nd_e_all[files_index] = nd_e_x
#         time[files_index] = files_index*scale*dt
            
# if(scale == 1 and len(sys.argv)>1 ):
#     nd_e_all.tofile("oscylacje_" + folder + ".csv", sep = ",")

# if(scale == 1):
#     nd_e_all = np.loadtxt("oscylacje_" + folder +".csv", delimiter=",").reshape(num_files, dimensions[0])
#     time = np.linspace(0,dt*(num_files-1), num_files)
# print(time)
# x = np.arange(0, dimensions[0], 1)*mesh.GetSpacing()[0]


# fig, ax = plt.subplots() 
# fig.set_size_inches(13,5)
# time *= 1e6
# pcolor = ax.pcolormesh(time, x, nd_e_all.transpose(), cmap='viridis', shading='auto')
# fig.colorbar(pcolor, label=r'$n_{srednie}$')  # Pasek kolorów
# ax.set_ylabel('$x$ [m]')
# ax.set_xlabel(r'$t$ [$\mu$s]')
# fig.savefig("obrazki/średnia gęstość w x2.pdf", dpi = 350)
# fig.savefig("obrazki/średnia gęstość w x2.png", dpi = 350)
# plt.subplots_adjust(left=0.15, right=0.95, top=0.92, bottom=0.15)
# # plt.show()

# how_much = slice(int(50/scale), len(time) - int(50/scale), 1) 
# place = 17
# # place = int(dimensions[0]/2)
# nd_mid_x = nd_e_all[how_much,place]
# fft_result = np.fft.fft(nd_mid_x)
# freqs = np.fft.fftfreq(len(nd_mid_x), d = dt)
# fig, ax = plt.subplots()

# ax.plot(time[how_much], nd_mid_x)
# ax.grid()
 

# fig, ax = plt.subplots()
# fig.set_size_inches(7,7)

# ax.scatter(freqs/(1e6)*2*np.pi, np.abs(fft_result))
# ax.set_xlabel(r"$\omega$ M[$\frac{\text{rad}}{\text{s}}$]")
# ax.set_ylabel("[-]")
# # ax.set_xlim((0,50))
# # ax.set_ylim((0, 0.8e13))
# ax.grid()
# plt.subplots_adjust(left=0.15, right=0.95, top=0.92, bottom=0.15)

# fig.savefig("obrazki/fft.pdf", dpi = 350)
# fig.savefig("obrazki/fft.png", dpi = 350)
# plt.show()
#     # print(f"Number of points: {mesh.GetNumberOfPoints()}")
#     # print(f"Number of cells: {mesh.GetNumberOfCells()}")
#     # print(f"mesh: {point_data}")

# ################################# oryginał: ################################3333

import numpy as np
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import sys

folder = "results_langmuir_2"

if(len(sys.argv) > 1): 
    scale = int(sys.argv[1])
else:
    scale = 1

plt.rcParams["font.size"] = 14

# Wczytaj plik .vtp
reader = vtk.vtkXMLImageDataReader()
num = str(f"{0:05}")
reader.SetFileName(r"../results_langmuir_2/fields_"+num+".vti")
reader.Update()
mesh = reader.GetOutput()
dimensions = mesh.GetDimensions()
point_data = mesh.GetPointData()

for i in range(point_data.GetNumberOfArrays()):
    if(point_data.GetArrayName(i) == "nd.e-"):
        nd_e_index = i

# num_files = 3004+1
num_files = int(1500/scale)
file_step = 30*scale
nd_e_all = np.zeros((num_files, dimensions[0]))
time = np.zeros(num_files)
dt = 20e-9
print(nd_e_all.shape)
print(time.shape)

if(len(sys.argv) > 1):
    for files_index in range(num_files):
        reader = vtk.vtkXMLImageDataReader()
        num = str(f"{files_index*file_step:05}")
        print(num)

        reader.SetFileName(r"../results_langmuir_2/fields_"+num+".vti")
        reader.Update()

        # Pobierz dane
        mesh = reader.GetOutput()
        # print(mesh)
        # print(f"mesh: {mesh}")
        # print(f"dimensions: {dimensions}")

        point_data = mesh.GetPointData()
        nd_e = point_data.GetArray(nd_e_index)
        # print(f"nd.e-: {nd_e}")
        nd_e_numpy = vtk_to_numpy(nd_e).reshape(dimensions)
        # print(nd_e_numpy)

        nd_e_x = np.zeros(dimensions[0])
        for k in range(dimensions[2]):
            nd_e_x[k] += nd_e_numpy[int(dimensions[0]/2)][int(dimensions[1]/2)][k]
        nd_e_x /= (dimensions[1]*dimensions[2])
        nd_e_all[files_index] = nd_e_x
        time[files_index] = files_index*scale*dt

if(scale == 1 and len(sys.argv)>1 ):
    nd_e_all.tofile("oscylacje_" + folder + ".csv", sep = ",")

if(scale == 1):
    nd_e_all = np.loadtxt("oscylacje_" + folder +".csv", delimiter=",").reshape(num_files, dimensions[0])
    time = np.linspace(0,dt*(num_files-1), num_files)
print(time)
x = np.arange(0, dimensions[0], 1)*mesh.GetSpacing()[0]


fig, ax = plt.subplots() 
fig.set_size_inches(13,5)
time *= 1e6
pcolor = ax.pcolormesh(time, x, nd_e_all.transpose(), cmap='viridis', shading='auto')
fig.colorbar(pcolor, label=r'$n_{srednie}$')  # Pasek kolorów
ax.set_ylabel('$x$ [m]')
ax.set_xlabel('$t$ [$\mu$s]')
fig.savefig("obrazki/średnia gęstość w x2.pdf", dpi = 350)
fig.savefig("obrazki/średnia gęstość w x2.png", dpi = 350)
plt.subplots_adjust(left=0.15, right=0.95, top=0.92, bottom=0.15)
# plt.show()

nd_mid_x = nd_e_all[:,int(3*len(x)/4)]
fft_result = np.fft.fft(nd_mid_x)
freqs = np.fft.fftfreq(len(nd_mid_x), d = dt)
fig, ax = plt.subplots()

ax.plot(time[:], nd_mid_x)
ax.grid()
 

fig, ax = plt.subplots()
fig.set_size_inches(7,7)

ax.scatter(freqs/(1e6)*2*np.pi, np.abs(fft_result))
ax.set_xlabel(r"$\omega$ M[$\frac{\text{rad}}{\text{s}}$]")
ax.set_ylabel("[-]")
ax.set_xlim((0,50))
ax.set_ylim((0, 0.8e13))
ax.grid()
plt.subplots_adjust(left=0.15, right=0.95, top=0.92, bottom=0.15)

fig.savefig("obrazki/fft.pdf", dpi = 350)
fig.savefig("obrazki/fft.png", dpi = 350)
plt.show()
    # print(f"Number of points: {mesh.GetNumberOfPoints()}")
    # print(f"Number of cells: {mesh.GetNumberOfCells()}")
    # print(f"mesh: {point_data}")
    
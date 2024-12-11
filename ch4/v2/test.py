# import sys
# import matplotlib.pyplot as plt
# import numpy as np
# from mpl_toolkits.mplot3d import Axes3D

# def parse_input(args):
#     # Convert command line arguments to arrays of floats
#     x1 = np.array([float(args[0]), float(args[1]), float(args[2])])
#     x2 = np.array([float(args[3]), float(args[4]), float(args[5])])
#     x3 = np.array([float(args[6]), float(args[7]), float(args[8])])
#     n = np.array([float(args[9]), float(args[10]), float(args[11])])
#     return x1, x2, x3, n

# # Parse command line input (excluding script name)
# x1, x2, x3, n = parse_input(sys.argv[1:])

# # Sphere parameters
# sphere_center = np.array([0, 0, 0.15])
# sphere_radius = 0.05

# # Create figure
# fig = plt.figure(figsize=(8, 8))
# ax = fig.add_subplot(111, projection='3d')
# #ax.set_box_aspect([1, 1, 1])  # Aspect ratio

# # Plot points and vectors
# ax.scatter(*x1, color='blue', label="x1", s=50)
# ax.scatter(*x2, color='green', label="x2", s=50)
# ax.scatter(*x3, color='red', label="x3", s=50)

# # Plot vector n originating from x3
# ax.quiver(*x3, *n, color='purple', length=0.1, normalize=True, label="Vector n")

# # Plot sphere with increased transparency
# u = np.linspace(0, 2 * np.pi, 100)
# v = np.linspace(0, np.pi, 50)
# x_sphere = sphere_center[0] + sphere_radius * np.outer(np.cos(u), np.sin(v))
# y_sphere = sphere_center[1] + sphere_radius * np.outer(np.sin(u), np.sin(v))
# z_sphere = sphere_center[2] + sphere_radius * np.outer(np.ones(np.size(u)), np.cos(v))
# ax.plot_surface(x_sphere, y_sphere, z_sphere, color='cyan', alpha=0.2, edgecolor='w')  # Increased transparency

# # Set labels and legend
# ax.set_xlabel("X-axis")
# ax.set_ylabel("Y-axis")
# ax.set_zlabel("Z-axis")
# ax.legend()

# # Set equal scaling for all axes
# max_range = np.array([x1, x2, x3, sphere_center + sphere_radius]).ptp(axis=0).max() / 2
# mid_x = (x1[0] + x2[0] + x3[0] + sphere_center[0]) / 4
# mid_y = (x1[1] + x2[1] + x3[1] + sphere_center[1]) / 4
# mid_z = (x1[2] + x2[2] + x3[2] + sphere_center[2]) / 4

# ax.set_xlim(mid_x - max_range, mid_x + max_range)
# ax.set_ylim(mid_y - max_range, mid_y + max_range)
# ax.set_zlim(mid_z - max_range, mid_z + max_range)

# plt.show()


import sys
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def parse_input(args):
    # Convert command line arguments to arrays of floats
    x1 = np.array([float(args[0]), float(args[1]), float(args[2])])
    x2 = np.array([float(args[3]), float(args[4]), float(args[5])])
    x3 = np.array([float(args[6]), float(args[7]), float(args[8])])
    n = np.array([float(args[9]), float(args[10]), float(args[11])])
    return x1, x2, x3, n

# Parse command line input (excluding script name)
x1, x2, x3, n = parse_input(sys.argv[1:])

# Cuboid parameters
cuboid_center = np.array([0, 0, 0])
cuboid_sides = np.array([0.1, 0.2, 0.3]) / 2  # Half-extents for each side

# Calculate cuboid vertices
vertices = np.array([
    [cuboid_center[0] - cuboid_sides[0], cuboid_center[1] - cuboid_sides[1], cuboid_center[2] - cuboid_sides[2]],
    [cuboid_center[0] + cuboid_sides[0], cuboid_center[1] - cuboid_sides[1], cuboid_center[2] - cuboid_sides[2]],
    [cuboid_center[0] + cuboid_sides[0], cuboid_center[1] + cuboid_sides[1], cuboid_center[2] - cuboid_sides[2]],
    [cuboid_center[0] - cuboid_sides[0], cuboid_center[1] + cuboid_sides[1], cuboid_center[2] - cuboid_sides[2]],
    [cuboid_center[0] - cuboid_sides[0], cuboid_center[1] - cuboid_sides[1], cuboid_center[2] + cuboid_sides[2]],
    [cuboid_center[0] + cuboid_sides[0], cuboid_center[1] - cuboid_sides[1], cuboid_center[2] + cuboid_sides[2]],
    [cuboid_center[0] + cuboid_sides[0], cuboid_center[1] + cuboid_sides[1], cuboid_center[2] + cuboid_sides[2]],
    [cuboid_center[0] - cuboid_sides[0], cuboid_center[1] + cuboid_sides[1], cuboid_center[2] + cuboid_sides[2]]
])

# Define cuboid faces using vertices
faces = [
    [vertices[0], vertices[1], vertices[2], vertices[3]],  # Bottom face
    [vertices[4], vertices[5], vertices[6], vertices[7]],  # Top face
    [vertices[0], vertices[1], vertices[5], vertices[4]],  # Front face
    [vertices[2], vertices[3], vertices[7], vertices[6]],  # Back face
    [vertices[1], vertices[2], vertices[6], vertices[5]],  # Right face
    [vertices[0], vertices[3], vertices[7], vertices[4]]   # Left face
]

# Create figure
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot points and vector
ax.scatter(*x1, color='blue', label="x1", s=50)
ax.scatter(*x2, color='green', label="x2", s=50)
ax.scatter(*x3, color='red', label="x3", s=50)
ax.quiver(*x3, *n, color='purple', length=0.1, normalize=True, label="Vector n")

# Plot cuboid
ax.add_collection3d(Poly3DCollection(faces, color='cyan', alpha=0.2, edgecolor='k'))

# Set labels and legend
ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")
ax.legend()

# Set equal scaling for all axes
max_range = np.array([x1, x2, x3, cuboid_center + cuboid_sides]).ptp(axis=0).max() / 2
mid_x = (x1[0] + x2[0] + x3[0] + cuboid_center[0]) / 4
mid_y = (x1[1] + x2[1] + x3[1] + cuboid_center[1]) / 4
mid_z = (x1[2] + x2[2] + x3[2] + cuboid_center[2]) / 4

ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

plt.show()

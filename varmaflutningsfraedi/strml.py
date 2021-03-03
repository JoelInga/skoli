import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d

r = np.logspace(0,1,200)
theta = np.linspace(0,np.pi*5,100)

N_r = len(r)
N_theta = len(theta)

# Polar to cartesian coordinates
theta_matrix, r_matrix = np.meshgrid(theta, r)
x = r_matrix * np.cos(theta_matrix)
y = r_matrix * np.sin(theta_matrix)

m = 5

# Stream function for a magnetic dipole
psi = -2*r_matrix*np.cos(theta_matrix)-theta_matrix-1/(4*np.pi)*np.log(r_matrix)

contour_levels = 5 * np.sin(np.linspace(0, np.pi/2,40))**2.

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.set_aspect('auto')

ax.plot_surface(x, y, psi, rstride=8, cstride=8, alpha=0.3)
ax.contour(x, y, psi, colors='black',levels=contour_levels)
plt.show()
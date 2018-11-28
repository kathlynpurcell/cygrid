import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


def func(x, y):
    return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2

grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]

# generate some data
points = np.random.rand(1000, 2)
values = func(points[:,0], points[:,1])

gridded = griddata(points, values, (grid_x, grid_y), method='cubic')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
ax1.imshow(func(grid_x, grid_y).T, extent=(0,1,0,1), origin='lower')
ax1.plot(points[:,0], points[:,1], 'k.', ms=1)
ax1.set_title('Original')
ax2.imshow(gridded.T, extent=(0,1,0,1), origin='lower')
ax2.set_title('Cubic interpolation')
plt.show()
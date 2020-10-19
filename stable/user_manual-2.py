import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


def func(x, y):
    return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2

grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]

fig, axes = plt.subplots(2, 2, figsize=(12, 12))
for ax, (sigma, nsize) in zip(axes.flatten(), [
        (0.01, 1000), (0.01, 100000), (0.1, 1000), (0.1, 100000)
        ]):

    points = np.random.rand(nsize, 2)
    noise = np.random.normal(0, sigma, nsize)

    values = func(points[:,0], points[:,1]) + noise
    gridded = griddata(points, values, (grid_x, grid_y), method='cubic')

    ax.imshow(gridded.T, extent=(0,1,0,1), origin='lower', vmin=-0.3, vmax=0.3)
    ax.set_title('nsize: {:d}, sigma: {:.2f}'.format(nsize, sigma))

plt.show()
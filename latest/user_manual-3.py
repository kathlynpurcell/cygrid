import numpy as np
import matplotlib.pyplot as plt
import cygrid


def func(x, y):
    return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2

grid_x, grid_y = np.meshgrid(
    np.linspace(0, 1, 100), np.linspace(0, 1, 100)
    )

fig, axes = plt.subplots(2, 2, figsize=(12, 12))
for ax, (sigma, nsize) in zip(axes.flatten(), [
        (0.01, 1000), (0.01, 100000), (0.1, 1000), (0.1, 100000)
        ]):

    points = np.random.rand(2, nsize)
    noise = np.random.normal(0, sigma, nsize)

    values = func(points[0], points[1]) + noise

    gridder = cygrid.SlGrid(grid_x.flatten(), grid_y.flatten(), 1)
    kernelsize_fwhm = 0.05
    kernelsize_sigma = kernelsize_fwhm / 2.355
    sphere_radius = 3. * kernelsize_sigma

    gridder.set_kernel(
        'gauss1d',
        (kernelsize_sigma,),
        sphere_radius,
        kernelsize_sigma / 2.
        )
    gridder.grid(points[0], points[1], values[:, None])
    gridded = gridder.get_datacube().squeeze().reshape(grid_x.shape)

    ax.imshow(gridded, extent=(0,1,0,1), origin='lower', vmin=-0.3, vmax=0.3)
    ax.set_title('nsize: {:d}, sigma: {:.2f}'.format(nsize, sigma))

plt.show()
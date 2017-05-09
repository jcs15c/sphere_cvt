import spherical_utils
import spherical_cvt

N = 24
M = 1000
generators = spherical_utils.uniform_sample(N)
"""
ax = spherical_utils.init_sphere()

spherical_utils.sphere_points(ax, generators)
spherical_cvt.plot_delaunay(ax,generators)
"""

bins, energy = spherical_cvt.bin_points(generators, M)
print(energy)


for k in range(50):
     spherical_cvt.cvt_step(generators, M)
     bins, energy = spherical_cvt.bin_points(generators, M)

print(energy)

ax = spherical_utils.init_sphere()

spherical_cvt.plot_voronoi(ax, generators)

spherical_utils.disp_sphere(ax)

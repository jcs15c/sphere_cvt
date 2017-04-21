import spherical_utils
import spherical_cvt

N = 100
M = 1000
generators = spherical_utils.uniform_sample(N)
print(generators[0])
"""
ax = spherical_utils.init_sphere()

spherical_utils.sphere_points(ax, generators)
spherical_cvt.plot_delaunay(ax,generators)
"""

bins, energy = spherical_cvt.bin_points(generators, M)
print(energy)

for k in range(5):
     spherical_cvt.cvt_step(generators, M)
     bins, energy = spherical_cvt.bin_points(generators, M)

print(energy)

ax = spherical_utils.init_sphere()

#spherical_utils.sphere_points(ax, generators)
spherical_cvt.plot_voronoi(ax, generators)
spherical_utils.disp_sphere(ax)
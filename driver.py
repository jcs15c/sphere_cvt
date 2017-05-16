import spherical_utils
import spherical_cvt
from projections import *

N = 100
M = 10000

plt.figure()
coast = read_coast_data()

generators = spherical_utils.uniform_sample(N)

bins, energy = spherical_cvt.bin_points(generators, M)

print(energy)

for k in range(10):
     print(k)
     spherical_cvt.cvt_step(generators, M)
     bins, energy = spherical_cvt.bin_points(generators, M)

print(energy)
"""
ax = spherical_utils.init_sphere()
spherical_cvt.plot_voronoi(ax, generators)
spherical_utils.disp_sphere(ax)
"""
plot_coast_map(coast, mercator)
spherical_cvt.plot_voronoi_map(generators, mercator, True)

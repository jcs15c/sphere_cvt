import spherical_utils
import spherical_cvt
from projections import *

N = 100
M = 10000

"""
gens = []
for i in range(10):
    gens.append(spherical_utils.get_cartesian(59.9343, 30 + i))
for i in range(10):
    gens.append(spherical_utils.get_cartesian(59. + i, 30 + i))
for i in range(10):
    gens.append(spherical_utils.get_cartesian(59. + i, 30))
for i in range(10):
    gens.append(spherical_utils.get_cartesian(69, 30 + i))
for i in range(10):
    gens.append(spherical_utils.get_cartesian(59. + i, 40))
generators = gens

for j in range(20):
    for i in range(5):
        generators.append(gens[i])
"""
generators = spherical_cvt.llarray_sample(N)
sample = spherical_cvt.prop_sample(M)

coast = read_coast_data()
plot_coast_map(coast, winkel)
spherical_cvt.plot_voronoi_map(generators, winkel, True)
plt.axis([-3, 3, -2, 2])
plt.savefig('test_0.png', dpi=500)

for k in range(30):
    if k:    
        plt.cla()
        coast = read_coast_data()
        plot_coast_map(coast, winkel)
        spherical_cvt.plot_voronoi_map(generators, winkel, True)
        plt.axis([-3, 3, -2, 2])
        plt.savefig('test_' + str(k) + '.png', dpi=500)
    
    sample = spherical_cvt.prop_sample(M)    
    spherical_cvt.cvt_step(generators, M, sample)

plt.cla()
coast = read_coast_data()
plot_coast_map(coast, winkel)
spherical_cvt.plot_voronoi_map(generators, winkel, True)
plt.axis([-3, 3, -2, 2])
plt.savefig('test_final.png', dpi=500)


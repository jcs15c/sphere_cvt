import spherical_utils
import spherical_cvt
from projections import *

N = 100
M = 100000

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
generators = spherical_cvt.llarray_sample(100)

coast = read_coast_data()
plot_coast_map(coast, winkel)
spherical_cvt.plot_voronoi_map(generators, winkel, True)
plt.axis([-3, 3, -2, 2])
plt.savefig('asia_init.png', dpi=500)

"""
plt.cla()
coast = read_coast_data()
plot_coast_map(coast, mercator)
spherical_cvt.plot_voronoi_map(generators, mercator, True)
plt.savefig('mercator_0.png', dpi=500)

plt.cla()
coast = read_coast_data()
plot_coast_map(coast, albers)
spherical_cvt.plot_voronoi_map(generators, albers, True)
plt.savefig('albers_0.png', dpi=500)
"""
bins, energy = spherical_cvt.bin_points(generators, M)
print(energy)

for k in range(30):
    if k:    
        plt.cla()
        coast = read_coast_data()
        plot_coast_map(coast, winkel)
        spherical_cvt.plot_voronoi_map(generators, winkel, True)
        plt.axis([-3, 3, -2, 2])
        plt.savefig('asia_' + str(k) + '.png', dpi=500)
        """
        plt.cla()
        coast = read_coast_data()
        plot_coast_map(coast, mercator)
        spherical_cvt.plot_voronoi_map(generators, mercator, True)
        plt.savefig('mercator_5.png', dpi=500)
        
        plt.cla()
        coast = read_coast_data()
        plot_coast_map(coast, albers)
        spherical_cvt.plot_voronoi_map(generators, albers, True)
        plt.savefig('albers_5.png', dpi=500)  
        """
    spherical_cvt.cvt_step(generators, M)
    bins, energy = spherical_cvt.bin_points(generators, M)
    print(k, energy)

print(energy)

plt.cla()
coast = read_coast_data()
plot_coast_map(coast, winkel)
spherical_cvt.plot_voronoi_map(generators, winkel, True)
plt.axis([-3, 3, -2, 2])
plt.savefig('asia_final.png', dpi=500)

"""
plt.cla()
coast = read_coast_data()
plot_coast_map(coast, mercator)
spherical_cvt.plot_voronoi_map(generators, mercator, True)
plt.savefig('mercator_30.png', dpi=500)

plt.cla()
coast = read_coast_data()
plot_coast_map(coast, albers)
spherical_cvt.plot_voronoi_map(generators, albers, True)
plt.savefig('albers_30.png', dpi=500)  
"""
import spherical_utils
import spherical_cvt
from projections import *

N = 100
M = 10000
"""
gens = []
for i in range(10):
    gens.append(spherical_utils.get_cartesian(42, -82 + i))
for i in range(10):
    gens.append(spherical_utils.get_cartesian(42. + i, -72 + i))
for i in range(10):
    gens.append(spherical_utils.get_cartesian(42. + i, -82))
for i in range(10):
    gens.append(spherical_utils.get_cartesian(52, -82 + i))
for i in range(10):
    gens.append(spherical_utils.get_cartesian(42. + i, -72))
generators = gens

for j in range(20):
    for i in range(5):
        generators.append(gens[i])
"""
t = 44

f = open("trial_" + str(t) + "/test_energy.txt",'w')

generators = spherical_utils.uniform_sample(N)
sample = spherical_utils.get_lat_long_array()

"""
coast = read_coast_data()
plot_coast_map(coast, winkel)
spherical_cvt.plot_voronoi_map(generators, winkel, True)
plt.axis([-3, 3, -2, 2])
plt.savefig("trial_" + str(t) + "/test_0.png", dpi=500)
"""
c = 0
lams = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
for j in lams:
    old_e = 1000000
    for k in range(1):
          
        plt.cla()
        coast = read_coast_data()
        plot_coast_map(coast, winkel)
        spherical_cvt.plot_voronoi_map(generators, winkel, True)
        plt.axis([-3, 3, -2, 2])
        plt.savefig('trial_' + str(t) + '/' + str(c) + '_test_' + str(k) + '_' + str(j) + '.png', dpi=500)
        
        #sample = spherical_cvt.llarray_sample(M)    
        new_e = spherical_cvt.cvt_step(generators, M, j, sample)
        f.write(str(j) + ' ' + str(k) + ' ' + str(new_e) + '\n')
        print(k, old_e - new_e)
        c += 1
        #           # break
        old_e = new_e
            
    
plt.cla()
coast = read_coast_data()
plot_coast_map(coast, winkel)

spherical_cvt.plot_voronoi_map(generators, winkel, True)
plt.axis([-3, 3, -2, 2])
plt.savefig('trial_' + str(t) + '/test_final.png', dpi=500)

f.close()


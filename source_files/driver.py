import spherical_utils
import spherical_cvt
import samples
import densities
from projections import *

N = 100
M = 10000

folder = "../output/test"

f = open(folder + "/test_energy.txt",'w')
f.write("iteration\t\tenergy\t\tstandard dev\t\tmean\n")

g = open(folder + "/document.txt",'w')
g.write("Proportional sample for generators CHOSEN WITH FULL UNLOG'D DENSITY, Uniform sample points, Full, LOG'D density")
g.close()

density = densities.get_den_array_log()
density2 = densities.get_den_array()
generators = samples.prop_sample(N, density2)
sample = samples.uniform_sample(M)

stdevs = []
means = []

for k in range(30):  
    if k%5 == 0:
        plt.cla()
        coast = read_coast_data()
        plot_coast_map(coast, lambert_cylindrical)
        spherical_cvt.plot_voronoi_map(generators, lambert_cylindrical, True)
        plt.axis([-3, 3, -2, 2])
        plt.savefig(folder + '/iteration_' + str(k) + '.png', dpi=250)
    
    e, stdev, mean = spherical_cvt.cvt_step(generators, k, folder, sample, density)
    stdevs.append(stdev)
    means.append(mean)

    f.write(str(k) + '\t\t' + str(e) + '\t\t' + str(stdev) +'\t\t' + str(mean) + '\n')
        
    
plt.cla()
coast = read_coast_data()
plot_coast_map(coast, lambert_cylindrical)

spherical_cvt.plot_voronoi_map(generators, lambert_cylindrical, True)
plt.axis([-3, 3, -2, 2])
plt.savefig(folder + '/iteration_30.png', dpi=500)

f.close()

plt.cla()
plt.plot(stdevs)
plt.savefig(folder + '/stdevs.png')

plt.cla()
plt.plot(means)
plt.savefig(folder + '/means.png')

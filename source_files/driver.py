import sys
import spherical_utils
import spherical_cvt
import sample_density
from projections import *

its = 50

pts = int(sys.argv[2])
trial = int(sys.argv[3])
sampler = str(sys.argv[1])

folder = sampler + '_' + str(pts)

f = open("../output/" + folder + "/trial_" + str(trial) + "/document.txt", 'w')
f.write("iteration\t\tenergy\t\tstd dev\t\tmean\n")

generators = sample_density.uniform_sample(50)[0]

if sampler == "lebedev":
    sample, density = sample_density.get_lebedev()
elif sampler == "fibonacci":
    sample, density = sample_density.get_fibonacci(pts)
elif sampler == "monte_carlo":
    sample, density = sample_density.get_monte_carlo(pts)
elif sampler == "helical":
    sample, density = sample_density.get_helical(pts)

stds = []
mean = []

for it in range(its + 1):
    plt.cla()
    coast = read_coast_data()
    plot_coast_map(coast, lambert_cylindrical)
    spherical_cvt.plot_voronoi_map(generators, lambert_cylindrical, True)
    plt.axis([-3, 3, -2, 2])
    energy, mv = spherical_cvt.cvt_step(generators, sample, density)#, folder, trial)
    if it%5 == 0:
        plt.savefig("../output/" + folder + "/trial_" + str(trial) + "/trial_" + str(it) + ".png")
        plt.cla()
        plt.plot(mv)
        plt.savefig("../output/" + folder + "/trial_" + str(trial) + "/mv_" + str(it) + ".png")
    f.write(str(it)+"\t\t"+str(energy)+"\t\t"+str(np.std(mv))+"\t\t"+str(np.mean(mv)) + '\n')
    mean.append(np.mean(mv))
    stds.append(np.std(mv))

plt.cla()
plt.plot(mean)
plt.savefig("../output/" + folder + "/trial_" + str(trial) + "/means.png")

plt.cla()
plt.plot(stds)
plt.savefig("../output/" + folder + "/trial_" + str(trial) + "/stds.png")

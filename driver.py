import spherical_utils
import spherical_cvt

N = 12
M = 1000
generators = spherical_utils.uniform_sample(N)

spherical_utils.plot_on_sphere(generators)


bins, energy = spherical_cvt.bin_points(generators, M)
print(energy)

for k in range(10):
	spherical_cvt.cvt_step(generators, M)
	bins, energy = spherical_cvt.bin_points(generators, M)
	print(energy)

spherical_utils.plot_on_sphere(generators)

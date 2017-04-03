import spherical_utils
import math
import numpy as np

def cvt_step(generators, M):
	
	bins, energy = bin_points(generators, M)
	
	for i in range(0, len(generators)):
		generators[i] = spherical_utils.project_to_unit_sphere(
										compute_center_of_mass(bins[i]))
	
def compute_center_of_mass(b):
	
	x1 = 0;
	y1 = 0;
	z1 = 0;
	
	for x in b:
		x1 += x[0]/len(b)
		y1 += x[1]/len(b)
		z1 += x[2]/len(b)
		
	return np.array([x1,y1,z1])
		
def bin_points(generators, M):
	
	# create M uniform test points
	x = spherical_utils.uniform_sample(M)
	
	bins = [ [] for x in range(len(generators))]
	energy = 0
	
	for xi in x:
		
		d = math.inf
		bi = -1
		
		for i in range(0,len(generators)):
			
			di = spherical_utils.distance_euclidean(xi,generators[i])
			if (di < d):
				d = di
				bi = i
		
		bins[bi].append(xi)
		energy += spherical_utils.distance_euclidean(xi, generators[bi])
		
	return bins, energy

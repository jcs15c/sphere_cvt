import numpy as np
import spherical_utils
import spherical_cvt
import densities
import random
from math import *

def uniform_sample(N):
	 """ Returns N uniformly random points on the surface of a sphere of
	 radius 1 centered at the origin.
  
	 Args
		N -- number of random points
		
	 Returns 
		  x -- a list of numpy arrays containing the cartesian coordinates 
			   of the points
	
	 """
	
	
	 x = list();
	 while len(x) < N:
		  x1 = random.gauss(0,1)
		  x2 = random.gauss(0,1)
		  x3 = random.gauss(0,1)
		  
		  r = sqrt(x1**2 + x2**2 + x3**2)
		  
		  nx = np.array([x1/r, x2/r, x3/r])
		  
		  x.append(nx)
		  
	 return x

def land_sample(N):
	llarray = get_land_array()
	generators = []
	
	for i in range(N):
		generators.append(llarray[np.random.randint(1, len(llarray)) - 1])
		
	return generators

def prop_sample(N, density):
	llarray = get_lat_long_array()
	prop_samp = []
	max_den = np.max(density)
	
	while len(prop_samp) < N:
		#Take random 0-1 real
		test_val = np.random.rand()
		#Take random non-zero sample point
		pt = llarray[np.random.randint(1, len(llarray)) - 1]
		#Convert to lat-long and get density at that point, divide by max
		lat, lon = spherical_utils.get_lat_long(pt)
		den = density[int(-(lat - 90)),int(lon + 180)]/max_den
		if den > test_val:
			prop_samp.append(pt)
	
	return prop_samp

### Array of all Lat-Long points ###
def get_lat_long_array():
	llarray = []
	for lon in range(-180, 181):
		for lat in range(-90, 91):
			llarray.append(spherical_utils.get_cartesian(lat, lon))
		
	return llarray

def get_land_array():
	pop_density = np.genfromtxt('../data/densities/population_full.csv', delimiter=',')
	land = []
	for lon in range(-180, 181):
		for lat in range(-90, 91):
			if pop_density[int(-(lat - 90)),int(lon + 180)] > 0:
				land.append(spherical_utils.get_cartesian(lat, lon))

	return land

def get_lebedev():
	file = open('../data/misc/lebedev_071.txt', 'r')
	pts = []

	for line in file:
	   temp = line.strip().split()
	   temp[0], temp[1] = float(temp[0]), float(temp[1])
	   pts.append(spherical_utils.get_cartesian(temp[0],temp[1]))

	file.close()
	return pts

def central_sample(lat, long):
	gens = []
	for i in range(10):
		for j in range(10):
			gens.append(spherical_utils.get_cartesian(lat + i, long+j))

	return gens
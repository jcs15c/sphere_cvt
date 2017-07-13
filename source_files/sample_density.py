import spherical_utils
import matplotlib.pyplot as plt
import scipy.interpolate
import numpy as np
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

		  y -- a list of uniform densities (1)
	
	 """
		
	 x = list();
	 y = []
	 while len(x) < N:
		  x1 = random.gauss(0,1)
		  x2 = random.gauss(0,1)
		  x3 = random.gauss(0,1)
		  
		  r = sqrt(x1**2 + x2**2 + x3**2)
		  
		  nx = np.array([x1/r, x2/r, x3/r])
		  
		  x.append(nx)
		  y.append(1)

	 return x, y

def get_lebedev():
	file = open('../data/lebedev_071.txt', 'r')
	#np.genfromtxt('../data/lebedev_densities.csv', delimiter=',')
	pts = []
	lats = []
	lons = []

	for line in file:
	   temp = line.strip().split()
	   lon, lat = float(temp[0]), float(temp[1])
	   lat -= 90
	   lats.append(lat)
	   lons.append(lon)
	   pts.append(spherical_utils.get_cartesian(lat, lon))

	file.close()
	return pts, get_densities(lats, lons)

def get_fibonacci(n):
	pts = []
	lats = []
	lons = []

	full_pts = sphere_fibonacci_grid_points(n)
	for i in full_pts:
		pts.append(i)
		lat, lon = spherical_utils.get_lat_long(i)
		lats.append(lat)
		lons.append(lon)

	return pts, get_densities(lats, lons)

def get_helical(n):
	pts = []
	h = np.zeros(n)
	lat = np.zeros(n)
	lon = np.zeros(n)

	for i in range(len(h)):
		h[i] = -1 + 2*(i)/(n-1.0)
		lat[i] = np.degrees(np.arccos(h[i]))

	for i in range(len(h)):
		if i == 0 or i == n-1:
			lon[i] = 0
			lon[n-1] = 0
		else:
			lon[i] = lon[i-1] + 3.6/np.sqrt(n)*1/np.sqrt(1 - h[i]**2)

	lon = np.degrees(lon)%360
	lat -= 90
	lon -= 180

	for i in range(len(lon)):
		pts.append(spherical_utils.get_cartesian(lat[i], lon[i]))

	return pts, get_densities(lat.tolist(), lon.tolist())

def get_monte_carlo(n):
	pts = uniform_sample(n)[0]
	lats = []
	lons = []

	for i in pts:
		lat, lon = spherical_utils.get_lat_long(i)
		lats.append(lat)
		lons.append(lon)

	return pts, get_densities(lats, lons)

def get_densities(lats, lons):
	cons = -407649103380480

	with np.load("../data/pop_density.npz") as data:
		ds = data['arr_0']
		dens = []
		for i in range(len(lats)):
			if lats[i] > 85 or lats[i] < -60:
				dens.append(cons)
			else:
				x, y = ll_to_grid(lats[i], lons[i])
				d = ds[x,y]
				dens.append(d)
	return dens

def ll_to_grid(lat,lon):
	y = int(120*(lon+180))
	x = int(120*(lat+60))
	if x:
		x -= 1
	if y:
		y -= 1
	return -x, y

def sphere_fibonacci_grid_points ( ng ):
#  Parameters:
#
#    Input, integer NG, the number of points.
#
#    Output, real XG(3,N), the grid points.
#
	phi = ( 1.0 + np.sqrt ( 5.0 ) ) / 2.0

	theta = np.zeros ( ng )
	sphi = np.zeros ( ng )
	cphi = np.zeros ( ng )

	for i in range ( 0, ng ):
		i2 = 2 * i - ( ng - 1 ) 
		theta[i] = 2.0 * np.pi * float ( i2 ) / phi
		sphi[i] = float ( i2 ) / float ( ng )
		cphi[i] = np.sqrt ( float ( ng + i2 ) * float ( ng - i2 ) ) / float ( ng )

	xg = np.zeros ( ( ng, 3 ) )

	for i in range ( 0, ng ) :
		xg[i,0] = cphi[i] * np.sin ( theta[i] )
		xg[i,1] = cphi[i] * np.cos ( theta[i] )
		xg[i,2] = sphi[i]

	return xg

def density_plot():
	cons = -407649103380480

	latitude = np.linspace(90, -90, 180)
	longitude = np.linspace(-180, 180, 360)

	[lons, lats] = np.meshgrid(longitude, latitude)

	shape = lats.shape

	dens = np.zeros(shape)
  
	with np.load("../data/densities/pop_density.npz") as data:
		ds = data['arr_0']
		for i in range(shape[0]):
			for j in range(shape[1]):
				if lats[i,j] > 85 or lats[i,j] < -60:
					dens[i,j] = 0
				else:
					x, y = ll_to_grid(lats[i,j], lons[i,j])
					d = ds[x,y]
					if d < 0:
						dens[i,j] = 0
					else:
						dens[i,j] = d
	
	dens = np.log10(dens + 10e-10) 

	plt.figure()
	plt.contourf(lons, lats, dens,30) 
	plt.title("raw data")
	plt.xlabel("longitude")
	plt.ylabel("latitude")
	plt.axis("equal")
	plt.colorbar()
	plt.show()
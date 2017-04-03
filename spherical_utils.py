import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
	for i in range(N):
		x1 = random.gauss(0,1)
		x2 = random.gauss(0,1)
		x3 = random.gauss(0,1)
		
		r = sqrt(x1**2 + x2**2 + x3**2)
		
		x.append(np.array([x1/r, x2/r, x3/r]))
		
	return x

########################################################################
########################################################################
def get_lat_long(x):
	""" Calculates the latitude and longitude of a cartesian point.
		
	This is the geocentric latitude and longitude, i.e. the polar and
	aximuth angles.
	
	Args
		x -- a point on the surface in cartesian coordinates

	Returns
		latitude -- the latitude of the point (polar angle)
		longitude -- the longitude of the point (azimuth angle)
	"""
	
	r = sqrt(x[0]**2 + x[1]**2 + x[2]**2)
	longitude = degrees(atan2(x[1],x[0]))
	latitude = degrees(acos(x[2]/r))
	
	#wrap latitude to [-90,90)
	latitude = latitude - 90
	
	return latitude, longitude

########################################################################
########################################################################
def distance_euclidean(x, y):
	""" Returns the Euclidean distance between two points
	
	Args
		x, y -- two numpy arrays representing points in R^3
		
	Returns
		d -- the Euclidean distance between x and y
		
	"""
	
	return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)
	
########################################################################
########################################################################
def project_to_unit_sphere(x):
	""" Projects a point inside the unit sphere to a point on the 
	surface.
	
	Args
		x -- a numpy array representing the cartesian coordinates of 
			a point inside the unit sphere
			
	Returns
		p -- the projection on x onto the surface of the unit sphere
		
	"""
	
	# compute spherical coordinates of x
	r = sqrt(x[0]**2 + x[1]**2 + x[2]**2)
	azimuth = atan2(x[1],x[0])
	polar = acos(x[2]/r)
	
	# compute cartesian coordinates of a point p with same azimuth and
	# polar angles as x, but with radius 1
	x1 = cos(azimuth)*sin(polar)
	y1 = sin(azimuth)*sin(polar)
	z1 = cos(polar)
	
	p = np.array([x1,y1,z1])
	
	return p
	
########################################################################
########################################################################
def plot_on_sphere(x):
	
	phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
	xx = np.sin(phi)*np.cos(theta)
	yy = np.sin(phi)*np.sin(theta)
	zz = np.cos(phi)
	
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	ax.plot_surface(
		xx, yy, zz,  rstride=1, cstride=1, color='c', alpha=0.3, linewidth=0)
	
	for xi in x:
		ax.scatter(xi[0],xi[1],xi[2],color="k",s=20)

	ax.set_xlim([-1,1])
	ax.set_ylim([-1,1])
	ax.set_zlim([-1,1])
	ax.set_aspect("equal")
	plt.tight_layout()
	plt.show()

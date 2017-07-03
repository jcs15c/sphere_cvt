import numpy as np
from scipy.special import sph_harm as harmonic
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def compute_coefficients(m, n, f, latitude, longitude):
    """
    Computes the coefficient for the spherical harmonic of order m and 
    degree n for a given data set defined over a sphere.
    
    Args:
        m [int]              
                -- harmonic order, |m| <= n
        n [int]              
                -- harmonic degree, n >= 0 
        f [2d numpy array]   
                -- data defined on the surface of the sphere
        latitude [2d numpy array] 
                -- array of latitude values in radians [0,pi]
        longitude [2d numpy array] 
                -- array of longitude values in radians [0,2pi]
        
    Returns:
        c [complex] 
                -- the coefficient (m,n) in the spherical harmonic
                        expansion of f
    """
    
    c = 0
    daz = np.pi/latitude.shape[0]
    dpol = 2*np.pi/longitude.shape[1]
    
    # compute integral over sphere surface
    # I'm using Reimann sums here, I think trapezoid would be better
    c = np.sum(f[1:,1:] * np.conj(harmonic(m, n, longitude[1:,1:], latitude[1:,1:]))*np.sin(latitude[1:,1:]))*daz*dpol
            
    return c

def compute_data_harmonic(m, n, c, latitude, longitude):
    """
    Computes the harmonic term (m,n) in a spherical harmonic expansion
    
    Args:
        m [int]              
                -- harmonic order, |m| <= n
        n [int]              
                -- harmonic degree, n >= 0 
        c [complex]   
                -- coefficient for the (m,n) term in expansion
        latitude [2d numpy array] 
                -- array of latitude values in radians [0,pi]
        longitude [2d numpy array] 
                -- array of longitude values in radians [0,2pi]
                
    Returns:
        f [2d numpy array]
                -- harmonic term evaluated at all lat/long locations
    """
    
    print("m = ",m,"n = ",n,"c = ",c)
    f = c*harmonic(m, n, longitude, latitude)
            
    return f

def compute_rms(c, n_max, data, latitude, longitude):
    """
    Computes the root-mean squared difference between pointwise data and
    a spherical harmonic expansion of degree n
    
    Args:
        c [list, numpy array of complex]
            -- spherical harmonic coefficients
        n_max [int]
            -- maximum degree of expansion
        data [2d numpy array]
            -- array containing actual data
        latitude [2d numpy array] 
            -- array of latitude values in radians [0,pi]
        longitude [2d numpy array] 
            -- array of longitude values in radians [0,2pi]   
             
    Returns:
        rms [float]
            -- root-mean squared difference between data and smoothed
                approximation   
    """
    
    data_smooth = np.zeros(data.shape)
    for n in range(n_max):
        for m in range(-n, n + 1):
            data_smooth += compute_data_harmonic(m, n, c[n][m + n], latitude, longitude).astype(float)
            
    rms = np.sqrt(np.sum((data - data_smooth)**2)/data.size)
    
    return rms

def get_den_array(numco = np.inf):
	if np.isinf(numco):
		return np.genfromtxt('../data/densities/population_full_adj.csv', delimiter=',')
	try:
		density = np.genfromtxt('../data/densities/pop_full_' + str(numco) + '.csv', delimiter=',')
	except Exception:
		return get_new_array(numco)
	return density

def get_den_array_log(numco = np.inf):
	if np.isinf(numco):
		return np.genfromtxt('../data/densities_log/population_log_adj.csv', delimiter=',')
	try:
		density = np.genfromtxt('../data/densities/pop_full_' + str(numco) + '.csv', delimiter=',')
	except Exception:
		return get_new_array(numco)
	return density

def get_new_array_log(numco):
	population_log = np.genfromtxt("../data/densities_log/population_log_adj.csv", delimiter=",")

	[n_lat, n_long] = population_log.shape
	latitude = np.linspace(np.pi,0,n_lat)
	longitude = np.linspace(0, 2*np.pi, n_long)

	[longitude, latitude] = np.meshgrid(longitude, latitude)

	# compute harmonic coefficients
	c = [] # start with empty list

	# set maximum degree
	N_MAX = numco
	for n in range(N_MAX + 1):
	    
	    c.append([]) 
	    # append a list of coefficients of degree n
	    for m in range(-n, n+1):
	        c[n].append(compute_coefficients(m, n, population_log, latitude, longitude))
	    
	    c[n] = np.asarray(c[n]) # convert list to numpy array, makes indexing easier
	    

	# compute spherical harmonic approximation to data
	population_smooth = np.zeros(population_log.shape)
	for n in range(N_MAX + 1):
	    for m in range(-n, n + 1):
	        population_smooth += compute_data_harmonic(m, n, c[n][m + n], latitude, longitude).astype(float)

	# shift approximation up to remove negative values
	min_value = np.min(population_smooth)
	population_smooth -= min_value
	
	np.savetxt("../data/densities_log/pop_log_" + str(numco) + ".csv", population_smooth, delimiter=",")

	return population_smooth

def get_new_array(numco):
	latitude = np.linspace(0, np.pi, 181)
	longitude = np.linspace(0, 2*np.pi, 361)

	[longitude, latitude] = np.meshgrid(longitude, latitude)    #Upside Down 

	data = np.genfromtxt('../data/densities/population_full_adj.csv', delimiter=',') 
	
	c = [] # start with empty list

	# set maximum degree
	N_MAX = numco

	for n in range(numco):
		
		c.append([]) 
		# append a list of coefficients of degree n
		for m in range(-n, n+1):
			c[n].append(compute_coefficients(m, n, data, latitude, longitude))
		
		c[n] = np.asarray(c[n]) # convert list to numpy array, makes indexing easier

	# compute spherical harmonic approximation to data
	data_smooth = np.zeros(data.shape)
	for n in range(N_MAX):
		for m in range(-n, n + 1):
			data_smooth += compute_data_harmonic(m, n, c[n][m + n], latitude, longitude).astype(float)

	# shift approximation up to remove negative values
	min_value = np.min(data_smooth)
	data_smooth -= min_value

	np.savetxt("../data/densities/pop_full_" + str(numco) + ".csv", data_smooth, delimiter=",")

	return data_smooth

import numpy as np
import matplotlib.pyplot as plt

### For Plotting Earth Map ###
def read_coast_data():
    file = open('world_coast.csv', 'r')
    coast = []
    
    for line in file:
       temp = line.strip().split(',')
       temp[0], temp[1] = float(temp[0]), float(temp[1])
       coast.append(temp)

    file.close()
    return coast

def plot_coast_map(coast, project):
    i = 0
    tempx = []
    tempy = []
    
    for i in range(len(coast)):
        if not np.isnan(coast[i][0]):
            coast[i][0], coast[i][1] = project(coast[i][1], coast[i][0])
            
    for i in range(len(coast)):
       if not np.isnan(coast[i][0]):
          tempx.append(coast[i][0])
          tempy.append(coast[i][1])
       else:
           plt.plot(tempx, tempy, 'r-')
           tempx = []
           tempy = []
       
### Projection Equations ###

def mercator(lat, long):
    return np.radians(long), np.log(np.tan(np.pi/4 + np.radians(lat)/2))
    
def lambert_cylindrical(lat, long):
    return np.radians(long), np.sin(np.radians(lat))

def albers(lat, long):
    lat = np.radians(lat)
    long = np.radians(long)
    
    sp1 = np.pi/9
    sp2 = 5*np.pi/18
    n = 0.5 * (np.sin(sp1) + np.sin(sp2))
    theta = n * long
    C = np.cos(sp1)**2 + 2 * n * np.sin(sp1)
    rho = 1 / n * np.sqrt(C - 2 * n * np.sin(lat))
    rho_null = 1 / n * np.sqrt(C)
    
    return rho * np.sin(theta), rho_null - rho * np.cos(theta)

def winkel(lat, long):
    lat = np.radians(lat)
    long = np.radians(long)
    
    a = np.arccos(np.cos(lat) * np.cos(long/2))
    sp = np.arccos(2/np.pi)
    
    x = 0.5 * (long * np.cos(sp) + ((2 * np.cos(lat) * np.sin(long/2))/usinc(a)))
    y = 0.5 * (lat + np.sin(lat)/usinc(a))
    
    return x, y
    
def aitoff(lat, long):
    lat = np.radians(lat)
    long = np.radians(long)
    
    a = np.arccos(np.cos(lat) * np.cos(long/2))
    
    x = (2 * np.cos(lat) * np.sin(long/2))/usinc(a)
    y = np.sin(lat)/usinc(a)
    
    return x, y

### Helper functions ###

def usinc(a):
    return np.sin(a) / a
    
def plot_outline(project):
    if project == mercator:
        plt.axis([-np.pi, np.pi, -np.pi/2, np.pi/2])
        plt.axvline(x = -np.pi)
        plt.axvline(x =  np.pi)
        
    if project == lambert_cylindrical:
        plt.axis([-np.pi, np.pi, -1, 1])
        plt.plot([-np.pi, np.pi, np.pi, -np.pi, -np.pi],
                 [1, 1, -1, -1, 1], 'b-')

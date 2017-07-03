import numpy as np
import matplotlib.pyplot as plt
import spherical_utils

def projection_point(pt, project):
    lat, long = spherical_utils.get_lat_long(pt)
    x, y = project(lat, long)
    plt.plot(x, y, 'b.')

def projection_line(u, v, project):
    #u & v are cartesian points
    r = 50
    
    t = np.linspace(0, 1, r)      
    w = v - u        
    x = u[0] + t*w[0]
    y = u[1] + t*w[1]
    z = u[2] + t*w[2]

    #Projects points on line onto unit sphere
    for i in range(r):
        pt = np.array([x[i],y[i],z[i]])
        new_pt = spherical_utils.project_to_unit_sphere(pt)
    
        x[i] = new_pt[0]
        y[i] = new_pt[1]
        z[i] = new_pt[2]

    #Converts points to lat/long
    lats = []
    long = []

    for i in range(r):
        lat_long = spherical_utils.get_lat_long([x[i],y[i],z[i]])
        lats.append(lat_long[0])
        long.append(lat_long[1])

    #Converts to xy
    cart_x = []
    cart_y = []
    for i in range(r):
        tx, ty = project(lats[i], long[i])
        cart_x.append(tx)
        cart_y.append(ty)
 
    d = np.pi/8
    #if np.linalg.norm(np.asarray([cart_x[0], cart_y[0]]) - np.asarray([cart_x[-1], cart_y[-1]])) < d:
    #    plt.plot(cart_x, cart_y, 'k-')
    #else:
    for i in range(r-1):
        if np.linalg.norm(np.asarray([cart_x[i], cart_y[i]]) - np.asarray([cart_x[i+1], cart_y[i+1]])) < d:
            plt.plot(cart_x[i:i+2], cart_y[i:i+2], 'k-')


### For Plotting Earth Map ###
def read_coast_data():
    file = open('../data/misc/world_coast.csv', 'r')
    coast = []
    
    for line in file:
       temp = line.strip().split(',')
       temp[0], temp[1] = float(temp[0]), float(temp[1])
       coast.append(temp)

    file.close()
    return coast

def plot_coast_map(coast, project):
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

def plot_coast_map_3d(ax, coast):
    tempx = []
    tempy = []
    tempz = []
                
    for i in range(len(coast)):
       if not np.isnan(coast[i][0]):
          pt = spherical_utils.get_cartesian(coast[i][1], coast[i][0])
          tempx.append(pt[0])
          tempy.append(pt[1])
          tempz.append(pt[2])
       else:
           #print(tempx)
           ax.plot(tempx, tempy, tempz, color='r')
           
           tempx = []
           tempy = []
           tempz = []

       
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

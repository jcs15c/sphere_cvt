import spherical_utils
import math
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import projections

pop_density = np.genfromtxt('population_full.csv', delimiter=',')
pop_density /= np.max(pop_density)
llarray = spherical_utils.get_lat_long_array()

def llarray_sample(N):
    generators = []
    
    for i in range(N):
        generators.append(llarray[np.random.randint(1, len(llarray)) - 1])
        
    return generators
    
def cvt_step(generators, M, lam, samples = np.nan):
     
     bins, energy = bin_points(generators, M, samples)
     
     for i in range(0, len(generators)):
          generators[i] = spherical_utils.project_to_unit_sphere(
                                                  compute_center_of_mass(bins[i], lam))
     return energy
     
def prop_sample(N):
    llarray = spherical_utils.get_lat_long_array()
    prop_samp = []
    max_den = np.max(pop_density)
    
    while len(prop_samp) < N:
        #Take random 0-1 real
        test_val = np.random.rand()
        #Take random non-zero sample point
        pt = llarray[np.random.randint(1, len(llarray)) - 1]
        #Convert to lat-long and get density at that point, divide by max
        lat, lon = spherical_utils.get_lat_long(pt)
        den = pop_density[int(-(lat - 90)),int(lon + 180)]/max_den
        if den > test_val:
            prop_samp.append(pt)
    
    return prop_samp
    
def compute_center_of_mass(b, lam):
     #if len(b) == 0:
     #    return prop_sample(1)[0]
     
     x1 = 0
     y1 = 0
     z1 = 0
     m  = 0
     
     for x in b:
          lat, lon = spherical_utils.get_lat_long(x)
          den = pop_density[int(-(lat - 90)),int(lon + 180)]
          if den < 0:
              den = 0
              
          den = 1 - lam + lam * den

              
          x1 += x[0] * den  
          y1 += x[1] * den
          z1 += x[2] * den
          m  += den
          
     return np.array([x1,y1,z1])/m
          
def bin_points(generators, M, samples = np.nan):
     
     # create M uniform test points
     if np.isnan(samples).all():
         x = spherical_utils.uniform_sample(M)
     else:
         x = samples
     
     bins = [ [] for x in range(len(generators))]
     energy = 0
     
     for xi in x:
          
          d = math.inf
          bi = -1
          
          for i in range(0,len(generators)):
               #lat, lon = spherical_utils.get_lat_long(xi)
               #den = pop_density[int(-(lat - 90)), int(lon + 180)] + 1e-10
               di = spherical_utils.distance_euclidean(xi,generators[i])# / den
               if (di < d):
                    d = di
                    bi = i
          
          bins[bi].append(xi)
          energy += spherical_utils.distance_euclidean(xi, generators[bi])
          
     return bins, energy

def plot_voronoi(ax, generators):
    tri_upp, tri_low, gen_upp, gen_low \
            = compute_delaunay(generators, True)
        
    pts_all, cen_all = recombine_del(gen_upp, gen_low, tri_upp, tri_low)
    
    ridges = calc_ridges(pts_all)
    
    for ridge in ridges:
        spherical_utils.sphere_line(ax, cen_all[ridge[0]], cen_all[ridge[1]], 'k')
    
def calc_ridges(pts):
    ridges = []
    
    for i in range(len(pts)):
        for j in range(len(pts)):
            if tri_adjacent(pts[i], pts[j]):
                ridges.append([i,j])

    return ridges
    
def tri_adjacent(tri1, tri2):
    found = 0    
    
    for pt in tri1:
        if (pt == tri2[0]).all():
            found += 1            
            break

    for pt in tri1:
        if (pt == tri2[1]).all():
            found += 1            
            break   

    for pt in tri1:
        if (pt == tri2[2]).all():
            found += 1   
            break

    if found == 2:
        return True
    return False        
    
def plot_delaunay(ax, generators):
    tri_upp, tri_low, gen_upp, gen_low \
            = compute_delaunay(generators, True)
        
    tri_all, cen_all = recombine_del(gen_upp, gen_low, tri_upp, tri_low)
    
    for tri in tri_all:
        spherical_utils.sphere_line(ax, tri[0], tri[1], 'k')
        spherical_utils.sphere_line(ax, tri[1], tri[2], 'k')
        spherical_utils.sphere_line(ax, tri[2], tri[0], 'k')
    
    spherical_utils.disp_sphere(ax)
    
def compute_delaunay(generators, full = 0):
    pts_upp = []
    pts_low = []
    gen_upp = []
    gen_low = []
    
    radius = np.pi - 0.1 
    
    for i in range(len(generators)):
        if np.arccos(np.dot([0,0,1],generators[i])) <= radius:
            pts_upp.append(spherical_utils.project_onto_upper_plane(generators[i]))
            gen_upp.append(generators[i])
        if np.arccos(np.dot([0,0,-1],generators[i])) <= radius:
            pts_low.append(spherical_utils.project_onto_lower_plane(generators[i]))
            gen_low.append(generators[i])

    tri_upp = Delaunay(pts_upp).simplices.copy()
    tri_low = Delaunay(pts_low).simplices.copy()
    
    #plot_2d_delaunay(pts_upp, tri_upp)
    
    tri_upp = remove_ill_tri(tri_upp, pts_upp, [0,0], radius)
    tri_low = remove_ill_tri(tri_low, pts_low, [0,0], radius)

    #plot_2d_delaunay(pts_upp, tri_upp)

    if full:
        return tri_upp, tri_low, gen_upp, gen_low
    return tri_upp, tri_low
    
def recombine_del(gen_upp, gen_low, tri_upp, tri_low):
    unique_tri = []    
    unique_cen = []

    for tri in tri_low:
        pt_tri = np.array([gen_low[tri[0]], gen_low[tri[1]], gen_low[tri[2]]])
        pt_cen = spherical_utils.cir_rad_center(pt_tri[0], pt_tri[1], pt_tri[2])[1]
        
        unique_tri.append(pt_tri)            
        unique_cen.append(pt_cen)            

    for tri in tri_upp:
        pt_tri = np.array([gen_upp[tri[0]], gen_upp[tri[1]], gen_upp[tri[2]]])
        pt_cen = spherical_utils.cir_rad_center(pt_tri[0], pt_tri[1], pt_tri[2])[1]

        found = False
        
        for i in range(len(unique_cen)):
            if (unique_cen[i] == pt_cen).all():
                found = True
                break
        if not found:
            unique_tri.append(pt_tri)            
            unique_cen.append(pt_cen)
         
    for i in range(len(unique_cen)):
        unique_cen[i] = spherical_utils.project_to_unit_sphere(unique_cen[i])
    
    return unique_tri, unique_cen

def plot_voronoi_map(generators, project, points = False):
    tri_upp, tri_low, gen_upp, gen_low \
            = compute_delaunay(generators, True)
        
    pts_all, cen_all = recombine_del(gen_upp, gen_low, tri_upp, tri_low)
    ridges = calc_ridges(pts_all)
    
    for ridge in ridges:
        projections.projection_line(cen_all[ridge[0]], cen_all[ridge[1]], project)

    if points:
        for generator in generators:
            projections.projection_point(generator, project)
    
    plt.axis('equal') 
    projections.plot_outline(project)
        
def plot_2d_delaunay(points, tri):
    plt.figure()
    points = np.asarray(points)
    plt.triplot(points[:,0], points[:,1], tri)
    plt.axis('equal')

    plt.ylim([-3,3])
    plt.xlim([-3,3])
    plt.show()
    
def remove_ill_tri(tri, pts, tcen, trad):
    i = len(tri)
    while i > 0:
        i += -1
        if not is_global(pts[tri[i,0]],pts[tri[i,1]],pts[tri[i,2]], tcen, trad):
            tri = np.delete(tri, i, axis=0)
    return tri
    
def is_global(p1, p2, p3, rcen, rrad):
    trad, tcen = spherical_utils.cir_rad_center(p1, p2, p3)
    #tcen = spherical_utils.project_to_unit_sphere(tcen)
    return ( (np.linalg.norm(rcen - tcen)) < rrad)
    
    
##############################################
    #####################################
    #####################################
##############################################
#Performs delaunay with only one tangent plane
def plot_delaunay_ill(ax, generators):
    pts = []    
    
    for i in range(len(generators)):
        pts.append(spherical_utils.project_onto_upper_plane(generators[i]))
        
    tri = Delaunay(pts)
    
    inds = tri.simplices.copy()
    
    for ind in inds:
        spherical_utils.sphere_line(ax, generators[ind[0]], generators[ind[1]])
        spherical_utils.sphere_line(ax, generators[ind[1]], generators[ind[2]])
        spherical_utils.sphere_line(ax, generators[ind[2]], generators[ind[0]])
        
    spherical_utils.disp_sphere(ax)

def plot_voronoi_tan(ax, generators, tangent):
    pts = []
    gen = []
    for i in range(len(generators)):
        if np.arccos(np.dot(tangent,generators[i])) <=  np.pi / 2:
            pts.append(spherical_utils.project_onto_tan_plane(generators[i], tangent))
            gen.append(generators[i])
   
    vor = Voronoi(pts)
    plot_internal_voronoi(vor)
    
    ridge_vertices = vor.ridge_vertices
    ridge_vertices = int_ridge(ridge_vertices)

    vertices_plane = vor.vertices
    vertices_sphere = np.zeros([vertices_plane.shape[0],3])
        
    for i in range(vertices_plane.shape[0]):
        vertices_sphere[i] = spherical_utils.project_onto_tan_sphere(vertices_plane[i], tangent)
       
    for ind in ridge_vertices:
        p1 = vertices_sphere[ind[0]]
        p2 = vertices_sphere[ind[1]]
        spherical_utils.sphere_line(ax,p1,p2)

def plot_internal_voronoi(vor):
    plt.figure()
    ridge_vertices = vor.ridge_vertices
    vertices = vor.vertices
    
    def int_ridge(ridge_vertices):
        i = len(ridge_vertices)
        while i > 0:
            i += -1
            if ridge_vertices[i][0] == -1 or \
               ridge_vertices[i][1] == -1:
                ridge_vertices = np.delete(ridge_vertices, i, axis=0)
        return ridge_vertices
    #print(ridge_vertices)
    ridge_vertices = int_ridge(ridge_vertices)
    #print(ridge_vertices)
    for ind in ridge_vertices:
        x = np.array([vertices[ind[0]][0],vertices[ind[1]][0]])
        y = np.array([vertices[ind[0]][1],vertices[ind[1]][1]])
        plt.plot(x, y, 'k-') 
        
def int_ridge(ridge_vertices):
    i = len(ridge_vertices)
    while i > 0:
        i += -1
        if ridge_vertices[i][0] == -1 or \
           ridge_vertices[i][1] == -1:
            ridge_vertices = np.delete(ridge_vertices, i, axis=0)
    return ridge_vertices



import spherical_utils
import math
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d


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

#Returns indicies of delaunay triangles
def compute_delaunay(generators, full = 0):
    pts_upp = []
    pts_low = []
    gen_upp = []
    gen_low = []
    
    for i in range(len(generators)):
        if np.arccos(np.dot([0,0,1],generators[i])) <= 3 * np.pi / 4:
            pts_upp.append(spherical_utils.project_onto_upper_plane(generators[i]))
            gen_upp.append(generators[i])
        if np.arccos(np.dot([0,0,-1],generators[i])) <= 3 * np.pi / 4:
            pts_low.append(spherical_utils.project_onto_lower_plane(generators[i]))
            gen_low.append(generators[i])

    tri_upp = Delaunay(pts_upp).simplices.copy()
    tri_low = Delaunay(pts_low).simplices.copy()
    
    #plot_2d_delaunay(pts_upp, tri_upp)
    
    tri_upp = remove_ill_tri(tri_upp, pts_upp, [0,0], 3 * np.pi / 4)
    tri_low = remove_ill_tri(tri_low, pts_low, [0,0], 3 * np.pi / 4)

    #plot_2d_delaunay(pts_upp, tri_upp)

    if full:
        return tri_upp, tri_low, gen_upp, gen_low
    return tri_upp, tri_low

def plot_voronoi(ax, generators):
    pts_upp = []
    
    for i in range(len(generators)):
        if np.arccos(np.dot([0,0,1],generators[i])) <= 3 * np.pi / 4:
            pts_upp.append(spherical_utils.project_onto_upper_plane(generators[i]))
       # if np.arccos(np.dot([0,0,-1],generators[i])) <= 3 * np.pi / 4:
    
    vor_upp = Voronoi(pts_upp)
    plot_internal_voronoi(vor_upp)
    
    rv_upp = vor_upp.ridge_vertices
    rv_upp = int_ridge(rv_upp)

    v_upp_plane = vor_upp.vertices
    v_upp_sphere = np.zeros([v_upp_plane.shape[0],3])
        
    for i in range(v_upp_plane.shape[0]):
        v_upp_sphere[i] = spherical_utils.project_onto_upper_sphere(v_upp_plane[i])
       
    for ind in rv_upp:
        p1 = v_upp_sphere[ind[0]]
        p2 = v_upp_sphere[ind[1]]
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
    
    
def plot_delaunay(ax, generators):
    tri_upp, tri_low, gen_upp, gen_low \
            = compute_delaunay(generators, True)
        
    for tri in tri_low:
        spherical_utils.sphere_line(ax, gen_low[tri[0]], gen_low[tri[1]], 'b')
        spherical_utils.sphere_line(ax, gen_low[tri[1]], gen_low[tri[2]], 'b')
        spherical_utils.sphere_line(ax, gen_low[tri[2]], gen_low[tri[0]], 'b')
    
    for tri in tri_upp:
        spherical_utils.sphere_line(ax, gen_upp[tri[0]], gen_upp[tri[1]], 'r')
        spherical_utils.sphere_line(ax, gen_upp[tri[1]], gen_upp[tri[2]], 'r')
        spherical_utils.sphere_line(ax, gen_upp[tri[2]], gen_upp[tri[0]], 'r')
    
    spherical_utils.disp_sphere(ax)

def plot_2d_delaunay(points, tri):
    plt.figure()
    points = np.asarray(points)
    plt.triplot(points[:,0], points[:,1], tri)
    plt.ylim([-5,5])
    plt.xlim([-5,5])
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





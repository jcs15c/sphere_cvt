# Spherical CVT by Population on Globe

This is a library to calculate and display Centroidal Voronoi Diagrams on the sphere, using global population data as a density function to create areas of roughly equal population. This is done using Python's Delaunay Triangulation functions and NASA's Socoeconomic Data and Application Center population data. 

## CVT on the sphere

While the calculations for computing the CVT clusters and centroids on a 3D surface are similar to the 2D case on a plane, or a 3D case in space, presenting the results in a comprehendable format requires additional functionality.

By utilizing a series of projections between the sphere and tangent planes at its poles, the vertices of the triangles of the Delaunay triangulation can be calculated and drawn on the sphere. Each point on the sphere is sorted and projected onto one of the two tangent planes, where the Triangulation is computed independently for the each half. The two halves are then recombined, and the corresponding delaunay ridges are drawn on the sphere.

THe Voronoi diagram can be computed from the circumcenters of each of the Delaunay triangles. Circumcenters for adjacent triangles are connected, and the resulting image is the Voronoi Diagram

```python
      generators = spherical_utils.uniform_sample(25)
      sample = spherical_utils.uniform_sample(1000)

      ax = spherical_utils.init_sphere()
      spherical_utils.sphere_points(ax, generators)
      spherical_utils.disp_sphere(ax)

      cx = spherical_utils.init_sphere()
      spherical_cvt.plot_voronoi(cx, generators)
      spherical_utils.disp_sphere(cx)

      dx = spherical_utils.init_sphere()
      spherical_cvt.plot_delaunay(dx, generators)
      spherical_utils.disp_sphere(dx)

      for i in range(5):
        #Third parameter governs density function, 0 for uniform
        spherical_cvt.cvt_step(generators, 1000, 0)

      bx = spherical_utils.init_sphere()    
      spherical_utils.sphere_points(bx, generators)
      spherical_utils.disp_sphere(bx)

      ex = spherical_utils.init_sphere()
      spherical_cvt.plot_voronoi(ex, generators)
      spherical_utils.disp_sphere(ex)

      fx = spherical_utils.init_sphere()
      spherical_cvt.plot_delaunay(fx, generators)
      spherical_utils.disp_sphere(fx) 
```
![Spherical CVT Example](https://github.com/jcs15c/sphere_cvt/blob/master/Spherical_CVT_Example.png "Spherical_CVT_Example")
 
## Removing Ill-Defined Voronoi Triangles
 
Because of the nature of the projection from a sphere to plane, some triangles are very narrow and are not true Delaunay triangles. These must be removed before the two haves are recombined, or else they will not overlap correctly. This is done by removing triangles whose circumcenters are outside of the radius of the circle of points assigned to each tangent plane. 

![Ill-Triange Example](https://github.com/jcs15c/sphere_cvt/blob/master/Ill_Triangle_Example.png "Ill_Triangle_Example")

## Coastal Map Data

Using existing MATLAB coast data, an outline of the world continents can be drawn on the 3D Sphere.

```python
      coast = read_coast_data()

      ax = spherical_utils.init_sphere()
      plot_coast_map_3d(ax, coast)
      spherical_utils.disp_sphere(ax)
```
![Spherical Coast Outline](https://github.com/jcs15c/sphere_cvt/blob/master/World_Map_3D.png "World_Map_3D")

# Spherical CVT by Population on Globe

This is a library to calculate and display Centroidal Voronoi Diagrams on the sphere, using global population data as a density function to create areas of roughly equal population. This is done using Python's Delaunay Triangulation functions and NASA's Socoeconomic Data and Application Center population data. 

## CVT on the sphere

While the calculations for computing the CVT clusters and centroids on a 3D surface are similar to the 2D case on a plane, or a 3D case in space, presenting the results in a comprehendable format requires additional functionality.


## Delaunay on the sphere

By utilizing a series of projections between the sphere and tangent planes at its poles, the vertices of the triangles of the Delaunay triangulation can be calculated and drawn on the sphere. Each point on the sphere is sorted and projected onto one of the two tangent planes, where the Triangulation is computed independently for the each half. 

```python
  #![]( "")
  
```

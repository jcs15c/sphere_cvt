from osgeo import gdal
import numpy as np

cons = -407649103380480

def ll_to_grid(lat,lon):
	y = int(120*(lon+180))
	x = int(120*(lat+60))
	if x:
		x -= 1
	if y:
		y -= 1
	return  17400 - x, y

xs = []
ys = []

ds = gdal.Open("gpw-v4-population-density_2015.tif").ReadAsArray()

file = open('lebedev_071.txt', 'r')
lebedev_densities = []

for line in file:
   	temp = line.strip().split()
   	lon, lat = float(temp[0]), float(temp[1])
   	lat -= 90

   	if lat > 85 or lat < -60:
   		lebedev_densities.append(cons)
   	else:
	   	x, y = ll_to_grid(lat,lon)
	   	d = ds[x,y]
		lebedev_densities.append(d)
		if d >= 0:
			xs.append(lon)
			ys.append(lat)


np.savetxt("lebedev_densities.csv", lebedev_densities, delimiter=",")
np.savetxt("xs.csv", xs, delimiter=",")
np.savetxt("ys.csv", ys, delimiter=",")


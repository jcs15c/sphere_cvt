import matplotlib.pyplot as plt
import numpy as np

xs = np.genfromtxt("xs.csv", delimiter=",")
ys = np.genfromtxt("ys.csv", delimiter=",")

plt.plot(xs,ys,'bo')
plt.axis([-180, 180, -90, 90])
plt.show()
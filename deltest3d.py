import numpy as np
from scipy.spatial import Delaunay
#from scipy.spatial import Voronoi, voronoi_plot_2d

points = np.random.rand(30, 3)
tri = Delaunay(points)
#or = Voronoi(points)
p = tri.points[tri.vertices]
n = tri.neighbors
for i in range(0,len(points)):
    ix=np.where(n==i) #indices of all triangles that contain the ith particle
    #print i, n[ix[0]] #print all the collections of triangles
    alltri=n[ix[0]]
    nn = np.setdiff1d(np.unique(np.concatenate(alltri)),np.array([i])) # all the neighboring indices
    nn = np.setdiff1d(nn,np.array([-1]))
    print i, len(nn),nn
#print n



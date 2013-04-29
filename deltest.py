import numpy as np
from scipy.spatial import Delaunay
#from scipy.spatial import Voronoi, voronoi_plot_2d

points = np.random.rand(30, 2)
tri = Delaunay(points)
#or = Voronoi(points)
p = tri.points[tri.vertices]
#print p
n = tri.vertices

for i in range(0,len(points)):
    ix=np.where(n==i) #indices of all triangles that contain the ith particle
    #print i, n[ix[0]] #print all the collections of triangles
    alltri=n[ix[0]]
    nn = np.setdiff1d(np.unique(np.concatenate(alltri)),np.array([i])) # all the neighboring indices
    print i, len(nn),nn

#print n


# Triangle vertices
A = p[:,0,:].T
B = p[:,1,:].T
C = p[:,2,:].T

# See http://en.wikipedia.org/wiki/Circumscribed_circle#Circumscribed_circles_of_triangles
# The following is just a direct transcription of the formula there
a = A - C
b = B - C

def dot2(u, v):
    return u[0]*v[0] + u[1]*v[1]

def cross2(u, v, w):
    """u x (v x w)"""
    return dot2(u, w)*v - dot2(u, v)*w

def ncross2(u, v):
    """|| u x v ||^2"""
    return sq2(u)*sq2(v) - dot2(u, v)**2

def sq2(u):
    return dot2(u, u)

cc = cross2(sq2(a) * b - sq2(b) * a, a, b) / (2*ncross2(a, b)) + C

# Grab the Voronoi edges
vc = cc[:,tri.neighbors]
vc[:,tri.neighbors == -1] = np.nan # edges at infinity, plotting those would need more work...


#print vc

lines = []
lines.extend(zip(cc.T, vc[:,:,0].T))
lines.extend(zip(cc.T, vc[:,:,1].T))
lines.extend(zip(cc.T, vc[:,:,2].T))

# Plot it
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

lines = LineCollection(lines, edgecolor='k')

plt.hold(1)
plt.plot(points[:,0], points[:,1], '.')
plt.plot(cc[0], cc[1], '*')
plt.gca().add_collection(lines)
plt.axis('equal')
plt.xlim(-0.1, 1.1)
plt.ylim(-0.1, 1.1)
plt.show()
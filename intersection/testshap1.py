import shapely
from shapely.geometry import Polygon
from pylab import plot,show
import numpy as np

a=np.array([[-1,1],[1,1],[1,-1],[-1,-1],[-1,1]])
p1=Polygon(a)

#p1 = Polygon([(0, 0), (1, 1), (1, 0)])
b=np.array([[0,2],[2,2],[2,0],[0,0],[0,2]])
p2 = Polygon(b)


plot(a[:,0],a[:,1])
plot(b[:,0],b[:,1])
x=p1.intersection(p2)
#print len(x)
show()
#x=p1.intersection(p2)

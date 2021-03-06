
import scipy
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay

import matplotlib.mlab as mlab

#inputfilelist=['uv_1.0.txt','uv_1.1.txt','uv_1.5.txt','uv_1.9.txt','uv_2.3.txt','uv_2.7.txt','uv_3.1.txt','uv_3.5.txt']
#aspectRatios=[1.,1.1,1.5,1.9,2.3,2.7,3.1,3.5]

aspectRatio=float(sys.argv[1])
inputfile='uv_'+str(aspectRatio)+'.txt'

#inputfilelist=['uv_3.9.txt']
#aspectRatios=[3.9]
N=50
L=512
R=L/2.
a=R
b=R
c=aspectRatio*R



#inputfile=sys.argv[1]
for i in range(0,1):
#for inputfile in inputfilelist:
    #inputfile=inputfilelist[i]
    f=open(inputfile,'r')
    lines=f.readlines()
    f.close
#lines=lines[len(lines)-51:len(lines)+1]
    u=[]
    v=[]
    un=np.zeros(N)
    vn=np.zeros(N)
    r=np.zeros(N)
    q=0
    for line in lines:
        these=line.split(',')
        u.append(float(these[0]))
        v.append(float(these[1]))
        un[q]=float(these[0])
        vn[q]=float(these[1])
        r[q]=float(these[3])
        q=q+1
    #print un
    x=a*np.cos(u)*np.cos(v)
    y=b*np.cos(u)*np.sin(v)
    z=c*np.sin(u)
    #print x
    #print y
    #print z
    #points=np.random.rand(30,3)
    #print points
    points=np.zeros((N,3))
    for i in range(0,N):
        points[i][0]=x[i]
        points[i][1]=y[i]
        points[i][2]=z[i]
    #print points
    tri=Delaunay(points)
    n=tri.vertices
    allDistances=np.array([])
    for i in range(0,len(points)):
        ix=np.where(n==i) #indices of all triangles that contain the ith particle
        #print i, n[ix[0]] #print all the collections of triangles
        alltri=n[ix[0]]
        #print alltri
        nn = np.setdiff1d(np.unique(np.concatenate(alltri)),np.array([i])) # all the neighboring indices
        nn = np.setdiff1d(nn,np.array([-1]))
        #print i, len(nn),nn
        #print "i=",i
        #print "ith point:",points[i]
        #print "neighbors:"
        #print points[nn]
        #print "differences:"
        #print (points[nn]-points[i])
        #print "distances:"
        r2n=(points[nn]-points[i])**2
        theseNN=np.array([])
        for t in r2n:
            d=np.sqrt(np.sum(t))
            #allDistances=np.append(allDistances,d)
            if d<(2*r[0]*1.5):
                theseNN=np.append(theseNN,d)
                
        allDistances=np.append(allDistances,len(theseNN))
        #print len(theseNN)
    
print allDistances
print r[0]*2
plt.hist(allDistances)
plt.show()
#print 2*r[0]
#plt.hist(allDistances,100)
#plt.show()
            #print dist
        #print "i=",i
        #ne=points[np.array(nn)]
        #print nn
        #for ne in nn:
        #    print ne



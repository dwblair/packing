
import scipy
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay

import matplotlib.mlab as mlab

#inputfilelist=['uv_1.0.txt','uv_1.1.txt','uv_1.5.txt','uv_1.9.txt','uv_2.3.txt','uv_2.7.txt','uv_3.1.txt','uv_3.5.txt']
#aspectRatios=[1.,1.1,1.5,1.9,2.3,2.7,3.1,3.5]

#aspectRatio=float(sys.argv[1])
#inputfile='uv_'+str(aspectRatio)+'.txt'

inputfile=str(sys.argv[1])
aspectRatio=1.
#inputfilelist=['uv_3.9.txt']
#aspectRatios=[3.9]
N=50
L=512
R=L/2.
a=R
b=R
c=aspectRatio*R

#inputfile='N_50_aspect_1.0_coords_000.txt'
f=open(inputfile,'r')
lines=f.readlines()
f.close()

maxT=1240 #arbitrary, needs to come in from last line in coords
t=0
l=0
timeStep=0
while timeStep<maxT:
    u=np.zeros(N)
    v=np.zeros(N) 
    r=np.zeros(N) 
    p=0
    points=np.zeros((N,2))
    for p in range(0,N):
        these=lines[l].split(',')
        u[p]=float(these[0])
        v[p]=float(these[1])
        r[p]=float(these[3])
        points[p][0]=float(these[0])
        points[p][1]=float(these[1])
        l=l+1
    these=lines[l].split(' ')
    timeStep=int(these[1])
    pressure=float(these[3])
    density=float(these[4])
    #print timeStep, density
    #print points
    tri=Delaunay(points)
    n=tri.vertices
    #print "timeStep=",timeStep
    psi6Ave=0.
    psi6Count=0
    for i in range(0,len(points)):
        ix=np.where(n==i)
        alltri=n[ix[0]]
        nn = np.setdiff1d(np.unique(np.concatenate(alltri)),np.array([i])) # all the neighboring indices
        nn = np.setdiff1d(nn,np.array([-1]))
        r2n=(points[nn]-points[i])**2
        theseNN=np.array([])
        #print nn
        finalNN=[]
        for j in range(0,len(nn)):
            d=np.sqrt(np.sum(r2n[j]))
            if d<.6:
                finalNN.append(int(nn[j]))
        #print i, finalNN
            #print nn[j],d
            #print nn[j]
        if len(finalNN)>2: #condition on having a psi6:
            neighbors=points[np.array(finalNN)]
            thisPoint=points[i]
            #print thisPoint,neighbors
            realPart=0.
            imPart=0.
            m=6.
            for k in range(0,len(neighbors)):
                dx=neighbors[k][0]-thisPoint[0]
                dy=neighbors[k][1]-thisPoint[1]
                r=np.sqrt(dx**2+dy**2)
                theta=np.arccos(dx/r)
                realPart=realPart+np.cos(theta*m)/len(neighbors)
                imPart=imPart+np.sin(theta*m)/len(neighbors)
            psiAbs2=realPart**2+imPart**2
            psi6=psiAbs2
            psi6Ave=psi6Ave+psi6
            psi6Count=psi6Count+1
    psi6Ave=psi6Ave/float(psi6Count)
    print timeStep, psi6Ave
            #print i, psi6

            #print points[np.array(finalNN)]
        #print i, theseNN
    l=l+1
plt.plot(u,v,'o')
plt.show()



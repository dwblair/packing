import numpy as np
import random
import time

def checkIntersections(i):
    neighborIndices=np.concatenate((np.arange(0,i,1),np.arange(i+1,N,1)),axis=0)
    #i.e. all the coords except i, for N**2
    #neighborIndices=np.arange(0,N)
    neighbors=coords[neighborIndices]
    dx=neighbors[:,0]-coords[i][0]
    dy=neighbors[:,1]-coords[i][1]
    r2=dx*dx+dy*dy
    intersectionIndices=r2<((2*r*rScale)**2) #i.e. less than a diameter
    numIntersections=len(r2[intersectionIndices])
    return numIntersections
    

def printOut():
    f=open('coords.txt','a')

    for i in range(0,N):
        thisline=str(coords[i][0])+","+str(coords[i][1])+","+str(radii[i]*rScale)+","+str(colors[i])+"\n"
        f.write(thisline)
            
    thisline="@ "+str(t)+" "+str(cellSize)+"\n"
    f.write(thisline)
    f.close()


timeGapForPrintout=1000

L=2**10 #system size, compatible with size of display in processing -- and must be a power of 2 for assignCell to work
N=4000 #number of particles
r=4. #particle radius

t=0
maxt=1000000

coords=np.zeros((N,2))+L/2.
colors=np.zeros(N)
radii=np.zeros(N)+r
stepSize=r*.1
rScale=1.
cellSize=10

# put coordinates in a grid
factor=1.8
x=factor*2*r
y=factor*2*r

for i in range(0,N):
    coords[i][0]=x;
    coords[i][1]=y;
    x=x+factor*2*r
    if x>(L-factor*2*r):
        x=factor*2*r
        y=y+factor*2*r
    radii[i]=r

timer_t=time.time()

# params for moving algorithm
mu = 0.
translationStep=1.

for t in range(0, maxt):

    #pick a random circle
    i=random.randint(0,N-1)
    #generate random displacements in x and y
    oldCoords=np.copy(coords[i])
    # generate random step
    dr=np.random.normal(mu,translationStep,2) #array of size 2, for each coordinate
    # move the particle
    coords[i]=coords[i]+dr

    if (checkIntersections(i))>0:
        coords[i]=oldCoords
    
    if t%timeGapForPrintout==0:
        printOut()
        print time.time()-timer_t
        timer_t=time.time()

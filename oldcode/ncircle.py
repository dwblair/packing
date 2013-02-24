import numpy as np
import random
import time


timeGapForPrintout=1000

L=512 #system size, compatible with size of display in processing -- and must be a power of 2 for assignCell to work
N=300 #number of particles
r=10 #particle radius

t=0
maxt=10000000

coords=np.zeros((N,2))+L/2.
radii=np.zeros(N)+r
stepSize=r*.1

# put coordinates in a grid
factor=1.2
x=factor*2*r
y=factor*2*r

for i in range(0,N):
    coords[i][0]=x;
    coords[i][1]=y;
    x=x+factor*2*r
    if x>(L-factor*2*r):
        x=factor*2*r
        radii[i]=r

timer_t=time.time()

for t in range(0, maxt):

    #pick a random circle
    i=random.randint(0,N-1)
    
    #move it
    coords[i][0]=coords[i][0]+random.normalvariate(0.,stepSize)
    coords[i][1]=coords[i][1]+random.normalvariate(0.,stepSize)

    if t%timeGapForPrintout==0:
        print time.time()-timer_t
        timer_t=time.time()

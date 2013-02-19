import numpy as np
import random
import time


class PackReplica:
    def __init__(self,L,N,timeGapForPrintout,outFileName,randomSeed):
        
        self.randomSeed=randomSeed
        self.L=L
        self.N=N
        self.outFileName=outFileName

        self.r=20.

        self.t=0
        self.maxt=1000000

        self.coords=np.zeros((N,2))+L/2.
        self.colors=np.zeros(N)
        self.cellNums=np.zeros(N)
        self.maxNumberParticlesPerCell=10  # assumes no more than 10 particles per cell! may need to modify if polydisperse ...
        #self.radii=np.zeros(N)+self.r
        #self.stepSize=self.r*.1
        self.rScale=1.
        self.cellSize=10
        self.MSD = np.zeros((N,2))

        # params for translation algorithm
        self.mu = 0.
        self.translationStep=10.
        self.expansionRatio=1.001
        self.ratioTranslateExpansion=100
        
        random.seed(self.randomSeed)


    def attemptTranslate(self,i): #random translate, particle i
        oldCoords=np.copy(self.coords[i])
        dr=np.random.normal(self.mu,self.translationStep,2) #array size 2 -- one random value for ea coord
        self.coords[i]=self.coords[i]+dr

        if (self.checkIntersections(i))>0:
            self.coords[i]=oldCoords

    def attemptExpansion(self):
        #particles=np.arange(0,self.N)
        #result=self.checkIntersections(particles)
        oldRScale=self.rScale
        self.rScale=self.rScale*self.expansionRatio
        intersections=0
        for i in range(0,self.N):
            intersections=intersections+self.checkIntersections(i)
            
        if intersections>0:
            self.rScale=oldRScale
        

    def sweep(self,numSweeps):
        #translation attempts
        for s in range(0,numSweeps):
            for p in range(0,self.ratioTranslateExpansion):
                i=random.randint(0,self.N-1)
                self.attemptTranslate(i)
            #expansion attempts
            self.attemptExpansion()
        self.t=self.t+1
        
    def checkIntersections(self,i):
        neighborIndices=np.concatenate((np.arange(0,i,1),np.arange(i+1,self.N,1)),axis=0)
        neighbors=self.coords[neighborIndices]
        dx=neighbors[:,0]-self.coords[i][0]
        dy=neighbors[:,1]-self.coords[i][1]
        r2=dx*dx+dy*dy
        # want to change this next to incorporate the diameter = radii[i]+radii[j], not global "r"
        #diams=self.radii*self.rScale
        
        intersectionIndices=r2<((2*self.r*self.rScale)**2) #i.e. less than a diameter
        numIntersections=len(r2[intersectionIndices])

        x=self.coords[i][0]
        y=self.coords[i][1]

        # in a box
        if x<(0+self.r*self.rScale) or x>(self.L-self.r*self.rScale):
            numIntersections=numIntersections+1
        if y<(0+self.r*self.rScale) or y>(self.L-self.r*self.rScale):
            numIntersections=numIntersections+1
        return numIntersections


        
    def initialize(self):

        #write out the params.in file
        p=open('params.in','w')
        thisline="#N,L,r\n"
        p.write(thisline)
        thisline=str(self.N)+","+str(self.L)+","+str(self.r)
        p.write(thisline)
        
        # write out a params.in file
        
        # place particles in a grid
        factor=1.5
        x=factor*2*self.r
        y=factor*2*self.r

        for i in range(0,self.N):
            self.coords[i][0]=x;
            self.coords[i][1]=y;
            x=x+factor*2*self.r
            if x>(self.L-factor*self.r):
                x=factor*2*self.r
                y=y+factor*2*self.r
            #self.radii[i]=self.r

        
    def printOut(self):
        f=open(self.outFileName,'a')

        for i in range(0,self.N):
            thisline=str(self.coords[i][0])+","+str(self.coords[i][1])+","+str(self.r*self.rScale)+","+str(self.cellNums[i])+"\n"
            f.write(thisline)

        thisline="@ "+str(self.t)+" "+str(self.cellSize)+"\n"
        f.write(thisline)
        f.close()


pr = PackReplica(L=2**9,N=10,timeGapForPrintout=10,outFileName='coords.txt',randomSeed=10)
pr.initialize()


for s in range(0,10000):
    pr.sweep(numSweeps=10)
    print s
    pr.printOut()

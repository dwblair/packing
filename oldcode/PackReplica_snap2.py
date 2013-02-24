import numpy as np
import random
import time


class PackReplica:
    def __init__(self,L,N,timeGapForPrintout,outFileName,pressure,randomSeed,initialDensity):

        
        self.pressure=pressure
        self.randomSeed=randomSeed
        self.L=L
        self.N=N
        self.outFileName=outFileName

        # assume particles are of radius "1"
        self.r=np.sqrt(initialDensity*L**2/(N*3.14159))
        self.V=self.N*np.pi/initialDensity

        self.density=initialDensity

        self.t=0
        self.maxt=1000000

        self.coords=np.zeros((N,2))+L/2.
        self.colors=np.zeros(N)
        self.cellNums=np.zeros(N)
        self.maxNumberParticlesPerCell=10  # assumes no more than 10 particles per cell! may need to modify if polydisperse ...
        #self.radii=np.zeros(N)+self.r
        #self.stepSize=self.r*.1
        #self.rScale=1.
        self.cellSize=10
        self.MSD = np.zeros((N,2))

        # params for translation algorithm
        self.mu = 0.
        self.translationStep=10.
        self.expansionStep=.001
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

        #print "attempt expand: r=",self.r

        
        dVScale=np.random.normal(self.mu,self.expansionStep)+1.
        dV=(dVScale-1.)*self.V

        
        
        expansionFlag=False
        if (dV < 0.):
            expansionFlag=True
        if (dV > 0.):
            #print dV, (1-np.exp(self.pressure*dV))
            #print dV, np.exp(-self.pressure*dV)
            p=random.random()
            boltz=min(1,np.exp(-self.pressure*dV))

            #print "pressure=",self.pressure,"; p=",p,"; boltz=",boltz 
            if p<boltz:
                expansionFlag=True
                #print "growing!"

        if expansionFlag==True:
            #print "Expansion!"
            oldVolume=self.V
            newVolume=self.V+dV
            
            #volumeRatio=newVolume/oldVolume
            self.V=newVolume
            oldR=self.r

            oldDensity=self.N*3.14159*self.r**2/self.L**2
            newDensity=self.N*3.14159/self.V #assumes particle radius 1
            #print "oldDensity=",oldDensity,"; newDensity=",newDensity

            self.r=np.sqrt(newDensity*self.L**2/(self.N*3.14159))
            #self.r=np.sqrt(self.L/self.V)
            #print "2=",np.sqrt(4.)
            #print "oldRatio=",np.sqrt(self.L/oldVolume)
            
            #print "oldR=",oldR,"; newR=",self.r
            
            #self.r=np.sqrt(oldR**2.)*oldVolume/newVolume
            
            intersections=0
            for i in range(0,self.N):
                intersections=intersections+self.checkIntersections(i)

            if intersections>0:
                self.r=oldR
                self.V=oldVolume
        
    def getDensity(self):
        self.density=self.N*3.14159*self.r**2/self.L**2
        return self.density
    
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
        
        intersectionIndices=r2<((2*self.r)**2) #i.e. less than a diameter
        numIntersections=len(r2[intersectionIndices])

        x=self.coords[i][0]
        y=self.coords[i][1]

        # in a box
        if x<(0+self.r) or x>(self.L-self.r):
            numIntersections=numIntersections+1
        if y<(0+self.r) or y>(self.L-self.r):
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
            thisline=str(self.coords[i][0])+","+str(self.coords[i][1])+","+str(self.r)+","+str(self.cellNums[i])+"\n"
            f.write(thisline)

        thisline="@ "+str(self.t)+" "+str(self.cellSize)+" "+str(self.getDensity())+"\n"
        f.write(thisline)
        f.close()


pr = PackReplica(L=2**9,N=10,timeGapForPrintout=10,outFileName='coords.txt',pressure=100.,randomSeed=10,initialDensity=.1)
pr.initialize()
print "L=",pr.L,"; r=",pr.r

for s in range(0,10000):
    pr.sweep(numSweeps=100)
    #print s, pr.r,pr.V/pr.L**2, pr.getDensity()
    print s, pr.r, pr.V
    pr.printOut()

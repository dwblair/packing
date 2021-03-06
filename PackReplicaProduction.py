import numpy as np
import random
import time


class PackReplica:
    def __init__(self,L,N,outFileNum,pressure,randomSeed,initialDensity):

        
        self.pressure=pressure
        self.randomSeed=randomSeed
        self.L=L
        self.N=N
        self.outFileNum=outFileNum

        # assume particles are of radius "1"
        self.r=np.sqrt(initialDensity*L**2/(N*3.14159))
        self.V=self.N*np.pi/initialDensity

        self.density=initialDensity

        self.t=0
        self.maxt=1000000

        self.coords=np.zeros((N,2))+L/2.

        # for MSD
        self.prevCoords=np.copy(self.coords)

        # for assessing equilibrium
        #self.prevDensity=self.density
        
        self.colors=np.zeros(N)
        self.cellNums=np.zeros(N)
        self.maxNumberParticlesPerCell=10  # assumes no more than 10 particles per cell! may need to modify if polydisperse ...
        #self.radii=np.zeros(N)+self.r
        #self.stepSize=self.r*.1
        #self.rScale=1.
        self.cellSize=10


        #variables for MSD calc
        self.prevCoords=np.zeros((N,2))

        #variables for assessing equilibrium
        self.equilibrated=False
        self.densitySTD=1. 

        # params for translation algorithm
        self.mu = 0.
        self.translationStep=10.
        self.expansionStep=.01
        
        self.ratioTranslateExpansion=100

        self.transAcceptRatio=0.
        self.expanAcceptRatio=0.
        self.expanRatioCount=0.
        self.transRatioCount=0.
        self.PIDfactor=0.1
        self.idealAcceptRatio=.4
        
        random.seed(self.randomSeed)


    def attemptTranslate(self,i): #random translate, particle i
        self.transRatioCount=self.transRatioCount+1 
        oldCoords=np.copy(self.coords[i])
        dr=np.random.normal(self.mu,self.translationStep,2) #array size 2 -- one random value for ea coord
        self.coords[i]=self.coords[i]+dr

        if (self.checkIntersections(i))>0:
            self.coords[i]=oldCoords
        else:
            self.transAcceptRatio=self.transAcceptRatio+1

    def attemptExpansion(self):
        #particles=np.arange(0,self.N)
        #result=self.checkIntersections(particles)

        #print "attempt expand: r=",self.r

        
        #dVScale=np.random.normal(self.mu,self.expansionStep)+1.
        
        dVScale=(2*random.random()-1)*self.expansionStep+1.
        #(2*random.random()-1.)*self.expansionStep
        dV=(dVScale-1.)*self.V
        #dV=dVScale*self.V-self.V

        
        
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
            self.expanRatioCount=self.expanRatioCount+1 
            oldVolume=self.V
            newVolume=self.V+dV
            
            #volumeRatio=newVolume/oldVolume
            self.V=newVolume
            oldR=self.r

            oldDensity=self.N*3.14159*self.r**2/self.L**2
            newDensity=self.N*3.14159/self.V #assumes particle radius 1
            #print "oldDensity=",oldDensity,"; newDensity=",newDensity

            if newDensity*self.L**2/(self.N*3.14159) <0:
                print "less than zero error!"

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
            else: # we accepted
                self.expanAcceptRatio=self.expanAcceptRatio+1
                #self.acceptRatioMeasureLength=1000.
                #self.acceptRatioCount=0
                
    
    def getDensity(self):
        self.density=self.N*3.14159*self.r**2/self.L**2
        return self.density
    
    def sweep(self,numSweeps):
        #translation attempts

        # assess density over range of numSweeps
        prevDensities=np.zeros(numSweeps)
        
        for s in range(0,numSweeps):
            for q in range(0,self.ratioTranslateExpansion): #makes sure we have a sweep of as many expansion attempts as we do the ratio between expansion attempts and translations
                for p in range(0,self.ratioTranslateExpansion):
                    i=random.randint(0,self.N-1)
                    self.attemptTranslate(i)
                #expansion attempts
                self.attemptExpansion()
            self.t=self.t+1 #however we define a sweep, we count t in 'sweeps'
            prevDensities[s]=self.getDensity()
            
        # update the step sizes using the PID loop
        self.updatePID() #modify the translation and expansion steps to get closer to ideal acceptance ratio

        # assess Equilibrium
        #print prevDensities
        self.densitySTD=np.std(prevDensities)
        #print self.t, self.PIDFactor, densitySTD
        
        
        
    def getMSD(self):
        # calcuates MSD, normalized by R**2
        #coordDiff=np.zeros((N,2))
        coordDiff=self.coords-self.prevCoords
        diffX=coordDiff[:,0]
        diffY=coordDiff[:,1]
        SD=diffX*diffX+diffY*diffY
        MSD=np.sum(SD)/float(self.N)/(self.r*self.r) #normalized by radius of particle
    
        #update prevCoords
        self.prevCoords=np.copy(self.coords)

        #return value
        return MSD
     
    def updatePID(self):
        
        self.transAcceptRatio=float(self.transAcceptRatio)/float(self.transRatioCount)
        self.expanAcceptRatio=float(self.expanAcceptRatio)/float(self.expanRatioCount)

        
        print self.t, self.PIDfactor, self.pressure, self.getDensity(), self.transAcceptRatio,self.expanAcceptRatio, self.translationStep, self.expansionStep, self.getMSD()

        #update the step sizes
        transDiff=self.transAcceptRatio-self.idealAcceptRatio
        expanDiff=self.expanAcceptRatio-self.idealAcceptRatio
        #self.PIDfactor=.1
        
        if transDiff>0: # accept ratio too high --> step too small
            self.translationStep=self.translationStep*(1.+self.PIDfactor)
            #print "accept too high; mod=",(1.+factor)
        if transDiff<0:  # accept ratio too low -->  step too big
            self.translationStep=self.translationStep*(1.-self.PIDfactor)

        if expanDiff>0: # accept ratio too high --> step too small
            self.expansionStep=self.expansionStep*(1.+self.PIDfactor)
        if  expanDiff<0: # accept ratio too low -->  step too big
            self.expansionStep=self.expansionStep*(1.-self.PIDfactor)
            
        #reset values
        self.transRatioCount=0
        self.expanRatioCount=0
        self.transAcceptRatio=0
        self.expanAcceptRatio=0
        
        
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

        #the coords file
        coordsOutFileName="coords"+str("%03d" % self.outFileNum)+".txt"
        densityOutFileName="density"+str("%03d" % self.outFileNum)+".txt"
        
        f=open(coordsOutFileName,'a')

        for i in range(0,self.N):
            thisline=str(self.coords[i][0])+","+str(self.coords[i][1])+","+str(self.r)+","+str(self.cellNums[i])+"\n"
            f.write(thisline)

        thisline="@ "+str(self.t)+" "+str(self.cellSize)+" "+str(self.pressure)+" "+str(self.getDensity())+"\n"
        f.write(thisline)
        f.close()

        f=open(densityOutFileName,'a')
        
        thisline=str(self.t)+" "+str(self.N)+" "+str(self.pressure)+" "+str(self.getDensity())+" "+str(self.densitySTD)+" "+str(self.getMSD())+" "+str(self.randomSeed)+" "+str(self.PIDfactor)+"\n"
        f.write(thisline)
        f.close()

        
#simNum=1
#pressure=.1
#pressureList=[10.,3.,2.,1.,.5,.1,0.1,.01,.001]
#NList=[5,6,7,8,10,11]
simNum=0
simNumMax=100
pressureList=[.01,.1, .2, .4, .8, 1., 5.,10.,12. ]
randomSeed=0
initialDensity=.1
N=11
L=2**9
tmax=2000
numSweeps=100
#while (pressure<10000):
#pr.initialize()

#for N in NList:
while simNum<simNumMax:
    
    pressure=pressureList[0]
    pr = PackReplica(L=L,N=N,outFileNum=simNum,pressure=pressure,randomSeed=randomSeed,initialDensity=initialDensity)
    pr.initialize()

    for pressure in pressureList:

        thisT=0

        pr.pressure=pressure

        #for s in range(0,10000)
        while (thisT<tmax):

            pr.sweep(numSweeps=numSweeps)
            #print pr.t, pr.N, pr.pressure, pr.PIDfactor, pr.getDensity(), pr.densitySTD, pr.getMSD()
            pr.printOut()

            thisT=thisT+numSweeps

        #pressure=pressure+.2
    randomSeed=randomSeed+1
    simNum=simNum+1

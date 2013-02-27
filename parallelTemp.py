import numpy as np
import random
import time
import os
import sys

import PackReplica


numSweeps=int(sys.argv[1])
tmax=int(sys.argv[2])
numSims=int(sys.argv[3])
N=int(sys.argv[4])
numReplicas=int(sys.argv[5])

simNum=0
invplist=np.linspace(40,.03,numReplicas)
plist=1./invplist
pressureList=plist[35:]

randomSeed=0
initialDensity=.05
L=2**9

# build up replicas
replicas=[]
for r in range(0,numReplicas):
    pressure=pressureList[0]
    simNum=r
    pr=PackReplica.PackReplica(L=L,N=N,outFileNum=simNum,pressure=pressure,randomSeed=randomSeed,initialDensity=initialDensity)
    pr.initialize()
    replicas.append(pr)

for p in pressureList:

    for r in range(0,numReplicas):
        


#for N in NList:
while simNum<numSims:
    
    pressure=pressureList[0]
    pr = PackReplica.PackReplica(L=L,N=N,outFileNum=simNum,pressure=pressure,randomSeed=randomSeed,initialDensity=initialDensity)
    pr.initialize()

    
    for pressure in pressureList:
        
        pr.pressure=pressure

        #equilibrate at first
        print "equilibrating ..."
        acceptRatiosGood=False

        # take a snapshot of the simulation before we let things evolve
        tSnap=pr.t
        pressureSnap=pr.pressure
        rSnap=pr.r
        VSnap=pr.V
        coordsSnap=np.copy(pr.coords)
        prevCoordsSnap=np.copy(pr.prevCoords)

             

        thisT=0
        while(acceptRatiosGood==False):
            pr.sweep(numSweeps=numSweeps)
            #pr.printOut()

            tr=pr.ratioStats[0] #translation accept ratio
            ex=pr.ratioStats[1] #expansion accept ratio
            print tr,ex, pr.MSD
            if (tr>.3 and tr<.5 and ex>.3 and ex<.5):
                acceptRatiosGood=True
                print "finished equilibrating in ", thisT, " steps"
                pr.requiredSweepsToAchieveOptimalAcceptRatio=thisT
                #thisT=tmax
            thisT=thisT+numSweeps

        #now begin actual simulation at this pressure value
        print "beginning serious sim!"

        #first retrieve original configurations
        pr.t=tSnap
        pr.pressure=pressureSnap
        pr.r=rSnap
        pr.V=VSnap
        pr.coords=np.copy(coordsSnap)
        pr.prevCoords=np.copy(prevCoordsSnap)

        #now carry out the sim
        thisT=0
        while (thisT<tmax):
            pr.sweep(numSweeps=numSweeps)
            pr.printOut()
            thisT=thisT+numSweeps
        
        #print pr.t, pr.PIDfactor, pr.pressure, pr.getDensity(), pr.transAcceptRatio,pr.expanAcceptRatio, pr.translationStep, pr.expansionStep, pr.getMSD(), pr.V

    randomSeed=randomSeed+1
    simNum=simNum+1

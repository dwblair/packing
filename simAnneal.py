import numpy as np
import random
import time
import os
import sys

import PackReplica



def copyReplicaSystem(r1,r2):
    
    r2.pressure=r1.pressure
    r2.L=r1.L
    r2.N=r1.N

    r2.r=r1.r
    r2.V=r1.V
    r2.density=r1.density
    r2.t=r1.t
    r2.prevCoords=np.copy(r1.prevCoords)

    r2.MSD=r1.MSD
    r2.mu=r1.mu

    r2.translationStep=r1.translationStep
    r2.expansionStep=r1.expansionStep
        
    r2.ratioTranslateExpansion=r1.ratioTranslateExpansion

    r2.transAcceptRatio=r1.transAcceptRatio
    r2.expanAcceptRatio=r1.expanAcceptRatio
    r2.expanRatioCount=r1.expanRatioCount
    r2.transRatioCount=r1.transRatioCount
    r2.PIDfactor=r1.PIDfactor
    r2.idealAcceptRatio=r1.idealAcceptRatio
    r2.ratioStats=r1.ratioStats
    r2.requiredSweepsToAchieveOptimalAcceptRatio=r1.requiredSweepsToAchieveOptimalAcceptRatio
        
    #density stats
    r2.highestDensity=r1.highestDensity

    #random seed
    #random.seed(self.randomSeed)






numSweeps=int(sys.argv[1])
tmax=int(sys.argv[2])
numSims=int(sys.argv[3])
N=int(sys.argv[4])
simNum=0
invplist=np.linspace(40,.03,50)
plist=1./invplist
pressureList=plist[35:]

randomSeed=0
initialDensity=.05

L=2**9

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

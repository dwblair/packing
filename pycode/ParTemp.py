
import numpy
from scipy import *
from pylab import *
import copy
import time

class ParTemp:
    def __init__(self, System, NumSystems, BetaMin, BetaMax, PrintInterval, StepsPerMove, NumMoves):
       self.System = System
       self.NumSystems = NumSystems
       self.NumMoves = NumMoves
       self.BetaMin = BetaMin
       self.BetaMax = BetaMax
       self.PrintInterval = PrintInterval
       self.StepsPerMove = StepsPerMove
       self.init_Beta()
       self.Systems = []
       for i in range(NumSystems):
           self.Systems.append(copy.deepcopy(System))
           self.Systems[i].Beta = self.Beta[i]
       return


    def init_Beta(self):
        norm = (1/self.BetaMin-1/self.BetaMax)/(self.NumSystems**2)
        self.Beta = zeros((self.NumSystems),dtype=float)
        for i in range(self.NumSystems):
            self.Beta[i] = 1/(norm * i*i + 1/self.BetaMax)
        return

    def move(self):
        #Emin = zeros((n),dtype=float)
        for i in range(self.NumMoves):
           for SysNum in range(self.NumSystems):
               System = self.Systems[SysNum]
               System.move(self.StepsPerMove)
               #Emin[SysNum] = System.getMinEnergy()

           RandNums = random((self.NumSystems))
           for SysNum in range(self.NumSystems-1):
               Beta1 = self.Beta[SysNum]
               Beta2 = self.Beta[SysNum+1]
               Energy1 = self.Systems[SysNum].Energy()
               Energy2 = self.Systems[SysNum+1].Energy()

               if RandNums[SysNum] < exp((Beta2-Beta1)*(Energy2-Energy1)):
                   self.Systems[SysNum], self.Systems[SysNum+1] = self.Systems[SysNum+1], self.Systems[SysNum] 

                   self.Systems[SysNum].Beta = Beta1
                   self.Systems[SysNum+1].Beta = Beta2

               if self.PrintInterval > 0 and i%self.PrintInterval == 0:
                   print "****iteration= ", i
                   for SysNum in range(self.NumSystems-1):
                       print "System ", System, "Beta = ", self.Beta[SysNum], " Energy = ", self.Systems[SysNum].Energy()
               if self.graph_on == 1:
                   #ion()
                   E = []
                   for i in range(self.NumSystems):
                       E.append(self.Systems[i].Energy())
                   cla()
                   plot(1/self.Beta,E)
                   show()
        return




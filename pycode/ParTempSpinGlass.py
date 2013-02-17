
import numpy
from scipy import *
from pylab import *
import time
from scipy import weave
from scipy.weave import converters
import ParTemp

class TwoDShortRangeSpinGlass:
    def __init__(self,Beta=1.0,size=32):
        self.Beta=Beta
        self.n = size
        self.s = ones((size,size),dtype=int)
        self.dim=2
        self.J = standard_normal((size,size,self.dim))
        
    def iterate(self):
       n = self.n
       r = random((n,n))
       s = self.s
       J = self.J
       beta = float(self.Beta)

       code = """
           #line 26 "ParTempSpinGlass.py"
           int mod = n-1;
           
           for(int x=0; x < n; x++) {
             for(int y=0; y < n; y++) {
                float field = J(x,y,0)*s(((x+1)&mod),y)+J(x,y,1)*s(x,((y+1)&mod))+J((x-1)&mod,y,0)*s(((x-1)&mod),y)+J(x,(y-1)&mod,1)*s(x,((y-1)&mod));
                float delta_e = s(x,y)*field;
                  if (exp(-beta*delta_e) > r(x,y)) {
                    s(x,y) = -s(x,y);
                    }
                 }
            }
               """
       weave.inline(code, ['J', 's', 'r', 'beta', 'n'], type_converters=converters.blitz, compiler = 'gcc')
       return

    def move(self,num_its):
        for i in range(num_its):
            self.iterate()
        return

    def Energy(self):
       s = self.s
       J = self.J
       n = self.n


       code = """
           #line 56 "ParTempSpinGlass.py"
           int x,y;
           int mod = n-1;
           float e_tot = 0;
           
           for(x=0; x < n; x++) {
             for(y=0; y < n; y++) {
                e_tot -= s(x,y)*(J(x,y,0)*s(((x+1)&mod),y)+J(x,y,1)*s(x,((y+1)&mod)));
                    }
                 }
            return_val=e_tot;
               """
       energy = weave.inline(code, ['J', 's', 'n'], type_converters=converters.blitz, compiler = 'gcc')
       return energy


sg = TwoDShortRangeSpinGlass(size=16)
partemp = ParTemp.ParTemp(System=sg, NumSystems=32, BetaMin=.3, BetaMax=10.0,PrintInterval=0,StepsPerMove=100, NumMoves=1000)
ion()
partemp.graph_on = 1;
partemp.move()
    



from numpy import *
import random

#### general brownian dynamics params #####

L=250 #system size, compatible with size of display in processing
N=20 #number of particles
r=20 #particle radius
t=0 #time
maxt=10000 #max time
stepSize=r/10. #random step size
dt=.1 #timestep
gamma=5. #friction coeff
kbt=2.5 # KbT
MAXFORCE=10000000. #in case something goes wrong
timeGapForPrintout=10 #number of timesteps to skip before printing out coordinates

##### lennard jones params #####
epsilon=100. #depth of lennard jones potential
s=r # 'radius' of potential well
s6=s**6 
s12=s6**2

#### gaussian random number parameters ######
sigma=1.
mu=0.

#coordinates array
coords=zeros((N,2))+L/2. #start all particles in the center

####### initially put coordinates in a grid #######
factor = 1.1
x=factor*2*r
y=factor*2*r

for i in range(0,N):
    coords[i][0]=x;
    coords[i][1]=y;
    x=x+factor*2*r
    if x>(L-factor*2*r):
        x=factor*2*r
        y=y+factor*2*r


############ run the simulation #########
for t in range(0,maxt): #time loop

    for i in range(0,N): #loop over all particles

        x=coords[i][0]
        y=coords[i][1]

        LJx=0.
        LJy=0.

        for j in range(0,N): #loop over all particles except the ith particle
            if j!=i:
                
                ##### get the distance between particles
                dx=coords[j][0]-coords[i][0]
                dy=coords[j][1]-coords[i][1]
                
                dr=sqrt(dx**2+dy**2)
                dr7=dr**7
                dr13=dr**13
               
                ###### calculate the LJ force in x and y

                ## note -- the neighbors need to be changed to reflect periodic boundaries (if not on compact surface)
                LJ=-24*epsilon*(2*s12/dr13 - s6/dr7)

                if (LJ>MAXFORCE):
                    LJ=MAXFORCE

                LJx=LJx+(dx/dr)*LJ
                LJy=LJy+(dy/dr)*LJ
         
        #### update the particle positions
        x=x+(dt/gamma)*LJx+sqrt(2*kbt*dt/gamma)*random.gauss(mu,sigma)
        y=y+(dt/gamma)*LJy+sqrt(2*kbt*dt/gamma)*random.gauss(mu,sigma)
        
        ##### enforce periodic boundaries (not on periodic surface)
        x=x%L
        y=y%L

        #put value back into coordinate array
        coords[i][0]=x
        coords[i][1]=y

        #output the particle positions
        if t%timeGapForPrintout==0:
            print coords[i][0],",",coords[i][1]

    #mark end of timestep
    if t%timeGapForPrintout==0:
        print "@", t

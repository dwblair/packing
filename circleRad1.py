from numpy import *
import random

#### general brownian dynamics params #####

L=250 #system size, compatible with size of display in processing
N=20 #number of particles
r=20 #particle radius
minR=1 #min particle radius
maxR=40 #max particle radius

t=0 #time
maxt=100 #max time
stepSize=r/10. #random step size
timeGapForPrintout=1 #number of timesteps to skip before printing out coordinates

#### gaussian random number parameters ######
sigma=1.
mu=0.

#coordinates array
coords=zeros((N,2))+L/2. #start all particles in the center

#radii array
radii=zeros(N)+r

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
    radii[i]=x*(maxR-minR)/L*1.+minR


############ run the simulation #########
for t in range(0,maxt): #time loop

    for i in range(0,N): #loop over all particles

        x=coords[i][0] 
        y=coords[i][1]


	xnew = x + (2*random.random()-1)*stepSize;
	ynew = y + (2*random.random()-1)*stepSize;


	# check overlaps
	overlaps=0

        for j in range(0,N): #loop over all particles except the ith particle
            if j!=i:
                
                ##### get the distance between particles
                #dx=coords[j][0]-coords[i][0]
                #dy=coords[j][1]-coords[i][1]
               
		dx = coords[j][0] - xnew;
		dy = coords[j][1] - ynew;
 
                dr=sqrt(dx**2+dy**2)
		if dr < (radii[i]+radii[j]):
			overlaps=overlaps+1
		
        	#inside a box
		if xnew > (L-radii[i]) or xnew < (0+radii[i]):
			overlaps=overlaps + 1
		if ynew > (L-radii[i]) or ynew < (0+radii[i]):
			overlaps=overlaps + 1 

        #### update the particle positions
	if overlaps==0:
		x=xnew
		y=ynew

        #x=x+(dt/gamma)*LJx+sqrt(2*kbt*dt/gamma)*random.gauss(mu,sigma)
        #y=y+(dt/gamma)*LJy+sqrt(2*kbt*dt/gamma)*random.gauss(mu,sigma)
        
        ##### enforce periodic boundaries (not on periodic surface)
        #x=x%L
        #y=y%L

        #put value back into coordinate array
        coords[i][0]=x
        coords[i][1]=y

        #output the particle positions
        if t%timeGapForPrintout==0:
            print coords[i][0],",",coords[i][1]

    #mark end of timestep
    if t%timeGapForPrintout==0:
        print "@", t

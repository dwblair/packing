from numpy import *
import random
#from __future__ import print_function

def checkOverlap(i,j):
    overlaps = 0
    #print "overlap for ",i," & ",j
    dx = coords[j][0] - coords[i][0]
    dy = coords[i][1] - coords[j][1]
    dr=sqrt(dx**2+dy**2)
    if dr < rScale*(radii[i]+radii[j]):
        overlaps=overlaps+1
    return overlaps

def checkOutsideBox(xnew,ynew):
    overlaps=0
    if xnew > (L-radii[i]) or xnew < (0+radii[i]):
        overlaps=overlaps + 1
    if ynew > (L-radii[i]) or ynew < (0+radii[i]):
        overlaps=overlaps + 1 
    return overlaps


#### general brownian dynamics params #####

L=500 #system size, compatible with size of display in processing
N=100 #number of particles
r=20 #particle radius
minR=1 #min particle radius
maxR=40 #max particle radius
BIAS=10

rScale=1.0
rStep=.001

t=0 #time
maxt=40000000 #max time
stepSize=.1*r #random step size
timeGapForPrintout=1000 #number of timesteps to skip before printing out coordinates

#### gaussian random number parameters ######
sigma=1.
mu=0.

#coordinates array
coords=zeros((N,2))+L/2. #start all particles in the center

#radii array
radii=zeros(N)+r

#neighbors array
neighbors=[]

####### initially put coordinates in a grid #######
factor = 1.1
x=factor*2*r
y=factor*2*r

# open the file
f=open('coords1.txt','w')

for i in range(0,N):
    coords[i][0]=x;
    coords[i][1]=y;
    x=x+factor*2*r
    if x>(L-factor*2*r):
        x=factor*2*r
        y=y+factor*2*r
    
    radii[i]=x*(maxR-minR)/L*1.+minR
    radii[i]=random.random()*10+10.
    theseNeighbors=arange(0,N,1)
    neighbors.append(theseNeighbors)

############ run the simulation #########
for t in range(0,maxt): #time loop


    overlaps=0


    # determine movetype
    moveType = random.randint(1,100)
    #moveType=1

    #print "moveType=",moveType

    if moveType != 1: # translation attempt
	
	i = random.randint(0,N-1)
	xold=coords[i][0]
	yold=coords[i][1]
        xnew = xold + (2*random.random()-1)*stepSize;
        ynew = yold + (2*random.random()-1)*stepSize;

	coords[i][0]=xnew
	coords[i][1]=ynew

        # check for being outside box
        overlaps=overlaps+checkOutsideBox(xnew,ynew)
        if (overlaps>0): #return particle to original position
            coords[i][0]=xold
            coords[i][1]=yold

        # if not outside box, check for overlaps with neighbors
	else: # check for overlaps in i's neighbors:
	    for n in neighbors[i]:
                if (i!=n):
                    overlaps=overlaps+checkOverlap(i,n)
                if (overlaps>0):
                    coords[i][0]=xold
                    coords[i][1]=yold
                    break  
		

    if moveType == 1: # grow / shrink attempt
        rScaleOld=rScale
        #BIASED -- used "BIAS*2" instead of "2":
	rScale = rScale + (2*BIAS*random.random()-1)*rStep
	#print "rScale=",rScale
	if (rScale > rScaleOld):
            for i in range(0,N):
                if overlaps>0:
                    rScale=rScaleOld
                    break
                for j in neighbors[i]:
                    if (i!=j):
                        overlaps=overlaps+checkOverlap(i,j)
                        if overlaps>0:
                            rScale=rScaleOld
                            break

    if t%timeGapForPrintout==0:
        for i in range(0,N):
            thisline=str(coords[i][0])+","+str(coords[i][1])+","+str(radii[i]*rScale)+"\n"
            f.write(thisline)
            
        thisline="@ "+str(t)+"\n"
	f.write(thisline)

f.close()



                        




        

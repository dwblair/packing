from numpy import *
import random
import time

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

def printOut(f):
    for i in range(0,N):
        thisline=str(coords[i][0])+","+str(coords[i][1])+","+str(radii[i]*rScale)+","+str(colors[i])+"\n"
        f.write(thisline)
            
    thisline="@ "+str(t)+" "+str(cellSize)+"\n"
    f.write(thisline)




def assignCells():
    #assumes that the system size L is an integer
    #takes the cell size to be the smallest power of 2 that is larger than or equal to the largest particle diameter in the system
    #this ensures that there will be an integer number of cell rows / columns in the system
    cellSize=2*ceil(max(radii)*rScale/2.)
    cellRows=L/cellSize
    power=ceil(log2(2*max(radii)*rScale))
    cellSize=2**power
    #print cellSize
    cellRows=L/cellSize
    cellContents=[]
    for i in range(0,N):
        cellContents.append([])

    for i in range(0,N):
        x=coords[i][0]
        y=coords[i][1]
        xCell=floor(x/cellSize)
        yCell=floor(y/cellSize)
        cellNum=yCell*cellRows+xCell
        cells[i]=cellNum
        colors[i]=cellNum
        thisCellContent=cellContents[int(cellNum)]
        thisCellContent.append(i)
        cellContents[int(cellNum)]=thisCellContent
        #print thisCellContent
        #cellContents[int(cellNum)]=cellContents[int(cellNum)].append(i)
    #print cellContents
    return cellSize

def assignCellContents():
    #assumes that the system size L is an integer
    #takes the cell size to be the smallest power of 2 that is larger than or equal to the largest particle diameter in the system
    #this ensures that there will be an integer number of cell rows / columns in the system
    cellSize=2*ceil(max(radii)*rScale/2.)
    cellRows=L/cellSize
    power=ceil(log2(2*max(radii)*rScale))
    cellSize=2**power
    #print cellSize
    cellRows=L/cellSize
    cellContents=[]
    for i in range(0,N):
        cellContents.append([])

    for i in range(0,N):
        x=coords[i][0]
        y=coords[i][1]
        xCell=floor(x/cellSize)
        yCell=floor(y/cellSize)
        cellNum=yCell*cellRows+xCell
        cells[i]=cellNum
        colors[i]=cellNum
        thisCellContent=cellContents[int(cellNum)]
        thisCellContent.append(i)
        cellContents[int(cellNum)]=thisCellContent
    #print cellContents 
    return cellContents 


def getNeighbors(i):
    neighbors=[]
    cellNum=int(cells[i])
    cellRows=int(L/cellSize)
    cellColumn=remainder(cellNum,cellRows)
    cellRow=floor(cellNum/cellRows)
   
    #its own cell
    thisContent=cellContents[cellNum]
    for n in thisContent:
        if n!=i:
            neighbors.append(n)
    #left
    if (cellColumn>0):
        thisContent=cellContents[cellNum-1]
        for n in thisContent:
            neighbors.append(n)
        #top left
        if (cellRow>0):
            thisContent=cellContents[cellNum-cellRows-1]
            for n in thisContent:
                neighbors.append(n)
        #bottom left
        if (cellRow<(cellRows-1)):
            thisContent=cellContents[cellNum+cellRows-1]
            for n in thisContent:
                neighbors.append(n)
    #right
    if (cellColumn<(cellRows-1)):
        thisContent=cellContents[cellNum+1]
        for n in thisContent:
            neighbors.append(n)
        #top right 
        if (cellRow>0):
            thisContent=cellContents[cellNum-cellRows+1]
            for n in thisContent:
                neighbors.append(n)
        #bottom right 
        if (cellRow<(cellRows-1)):
            thisContent=cellContents[cellNum+cellRows+1]
            for n in thisContent:
                neighbors.append(n)
        
    #top
    if (cellRow>0):
        thisContent=cellContents[cellNum-cellRows]
        for n in thisContent:
            neighbors.append(n)
    #bottom
    if (cellRow<(cellRows-1)):
        thisContent=cellContents[cellNum+cellRows]
        for n in thisContent:
            neighbors.append(n)
    return neighbors




#### general brownian dynamics params #####

L=512 #system size, compatible with size of display in processing -- and must be a power of 2 for assignCell to work
N=300 #number of particles
r=10 #particle radius
minR=1 #min particle radius
maxR=r #max particle radius
BIAS=2
cellSize=2

rScale=1.0
rStep=.001

t=0 #time
maxt=2000000 #max time
stepSize=.1*r #random step size
timeGapForPrintout=100 #number of timesteps to skip before printing out coordinates

#### gaussian random number parameters ######
sigma=1.
mu=0.

#coordinates array
coords=zeros((N,2))+L/2. #start all particles in the center

#radii array
radii=zeros(N)+r

#neighbors array
neighbors=[]

#cells array
cells=zeros(N)
cellContents=[]
for i in range(0,N):
    cellContents.append([])
#thisCell=cellContents[10]
#cellContents[3]=thisCell.append([4,2])
#print cellContents


cellNeighbors=[]


#colors array
colors=zeros(N)

####### initially put coordinates in a grid #######
factor = 1.2 
x=factor*2*r
y=factor*2*r

# open the file
f=open('coords.txt','w')

for i in range(0,N):
    coords[i][0]=x;
    coords[i][1]=y;
    x=x+factor*2*r
    if x>(L-factor*2*r):
        x=factor*2*r
        y=y+factor*2*r
    #radii[i]=x*(maxR-minR)/L*1.+minR
    #radii[i]=random.random()*10+10.
    radii[i]=r
    theseNeighbors=arange(0,N,1)
    neighbors.append(theseNeighbors)

timert=time.time()

############ run the simulation #########
for t in range(0,maxt): #time loop

    

    #print t
    #assign cells
    cellSize=assignCells()
    cellContents=assignCellContents()

    #print "cellRows",L/cellSize
    #print "cellContents:"
    #for i in range (0,N):
    #    print i, cellContents[i]
    #print "neighbors of 3",getNeighbors(3)
    #print cellContents

    #for i in range(0,N):
    #    print i, getNeighbors(i)

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
            for n in getNeighbors(i):
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
        if (rScale > rScaleOld):
            for i in range(0,N):
                if overlaps>0:
                    rScale=rScaleOld
                    break
                for j in getNeighbors(i):
                    if (i!=j):
                        overlaps=overlaps+checkOverlap(i,j)
                        if overlaps>0:
                            rScale=rScaleOld
                            break

    if t%timeGapForPrintout==0:
        printOut(f)
        print time.time()-timert
        timert=time.time()
        #for i in range(0,N):
        #    thisline=str(coords[i][0])+","+str(coords[i][1])+","+str(radii[i]*rScale)+"\n"
        #    f.write(thisline)
            
        #thisline="@ "+str(t)+"\n"

f.close()



                        




        

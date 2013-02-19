import sys
sys.path.append('/u/home/campus/dwblair/lib/python')
import Polygon, Polygon.IO, Polygon.Shapes
from numpy import *
import math
import random
import time
import cPickle as pickle


class Replica:

	def __init__(self):

		#print "initializing!"

		self.pressure=0.
		self.density=0.
		self.N=0
		
		self.L=0.  # set this; the rect_side is what gets modified
	
		self.replicaLabel=0

		self.volume=0.
 
		self.positions=[]
		#self.prev_positions=[]
		self.velocities=[]

		self.polygons=[]
		self.polygonPositions=[]
		self.polygonContoursInitial=[]
               
		self.polygonsJ=[]
		self.polygonPositionsJ=[]
                
		self.shapeTypes=[]

		self.translate_step=1.
		self.rotate_step=1.
		self.scale_step=.1

		self.p_translate=.4
		self.p_rotate=self.p_translate+.4

		self.translate_accept=0.
		self.translate_attempt=0

		self.rotate_accept=0.
		self.rotate_attempt=0

		self.scale_accept=0.
		self.scale_attempt=0

		self.accept_goal=.4

		self.K_c=.1

		self.stat_collect=100 # interval for collecting statistics	

		self.t=0
		self.sweeps=0.
		self.dVscale=0.
		self.dVabsolute=0.
		
		self.boundingBoxesInitial=[]

		#self.cellL=0.
		
		self.r2cutoff=0.
		self.r2neighborlist=[]
		self.addedR2Buffer=1.
		self.r2neighborUpdateCount=100
		self.r2cutoffUpdateCount=100

		self.ave_velocity=0.
		self.total_displacement=0.
		
		#self.initialSize=0.
		#self.finalSize=0.
		self.initialContours=[]


        def scalePolygons(self,scaleRatio):

            for pi in range(len(self.polygons)):
                (cx,cy)=self.polygons[pi].center()
                self.polygons[pi].scale(scaleRatio,scaleRatio,cx,cy)
                (cx,cy)=self.polygonsJ[pi].center()
                self.polygonsJ[pi].scale(scaleRatio,scaleRatio,cx,cy)

	
                
	def movePolygon(self,pi,x,y,a):
		
		pix=self.polygonPositions[pi][0]
		piy=self.polygonPositions[pi][1]
		pia=self.polygonPositions[pi][2]

                
		#print "x,y,a=",x,y,a


		da=a-pia
		dy=y-piy
		dx=x-pix

                #print "dx,dy,da=",dx,dy,da

                #cx=self.polygonCenters[pi][0]
                #cy=self.polygonCenters[pi][1]
                
                (cx,cy)=self.polygons[pi].center()


		#move the polygon
		self.polygons[pi].shift(dx,dy)

                (cx,cy)=self.polygons[pi].center()
		self.polygons[pi].rotate(da,cx,cy)
		
		#update polygonPositions
		newPosition=(x,y,a)
		self.polygonPositions[pi]=newPosition
		
	
	def movePolygonJ(self,pi,x,y,a):
		
		pix=self.polygonPositionsJ[pi][0]
		piy=self.polygonPositionsJ[pi][1]
		pia=self.polygonPositionsJ[pi][2]
		
		da=a-pia
		dy=y-piy
		dx=x-pix

                #cx=self.polygonCentersJ[pi][0]
                #cy=self.polygonCentersJ[pi][1]

                #cx=self.polygonsJ[pi].cx
                #cy=self.polygonsJ[pi].cy
                self.polygonsJ[pi].shift(dx,dy)
                (cx,cy)=self.polygonsJ[pi].center()

		#move the polygon
	
		self.polygonsJ[pi].rotate(da,cx,cy)
		
		#update polygonPositions
		newPosition=(x,y,a)
		self.polygonPositionsJ[pi]=newPosition

                
	def getTotalParticleArea(self):

		#print "getTotalParticleArea()"

		totalParticleArea=0.
		for i in range(self.N):
			thisShape=int(self.shapeTypes[i])
			totalParticleArea=totalParticleArea+self.polygons[thisShape].area()
		return totalParticleArea

	def updateNeighborList(self,addedBuffer):

		self.r2neighborlist=[]
		#theseNeighbors=[]
		#print "r2cutoff=",self.r2cutoff
		#print "updating neighbors"

		for j in range(self.N):
			theseNeighbors=[]
			for i in range(self.N):
				if i!=j:
					xi=self.positions[i][0]
					yi=self.positions[i][1]
					#ai=self.positions[i][2]

					xj=self.positions[j][0]
					yj=self.positions[j][1]
					#aj=self.positions[j][2]
					
					# periodic
					if (xi-xj) > (self.L/2.):
						xj=xj+self.L
					if (xi-xj) < (-self.L/2.):
						xj=xj-self.L
					if (yi-yj) > (self.L/2.):
						yj=yj+self.L
					if (yi-yj) < (-self.L/2.):
						yj=yj-self.L

					dx=xi-xj
					dy=yi=yj

					r2=dx**2+dy**2
					#print "r2=",r2
					
					if r2<(self.r2cutoff+addedBuffer):
						theseNeighbors.append(i)
			self.r2neighborlist.append(theseNeighbors)
						

        def findAnyOverlapsWithNeighbors(self,i):

                anyOverlap=0

		"""
		for j in self.r2neighborlist[i]:
			doesOverlap=self.testOverlapByIndex(i,j)
			if (doesOverlap==True):
				anyOverlap=anyOverlap+1

		"""
		for j in range(self.N):
			if i!=j:
				doesOverlap=self.testOverlapByIndex(i,j)
				if (doesOverlap==True):
					anyOverlap=anyOverlap+1
					#print "YES!"
					
		

		if (anyOverlap>0):
			return True
		else:
			return False


	def randMove(self):
		for i in range(self.N):
			maxL=self.L/2
			x=(2*random.random()-1.)*maxL
			#if i==0:
				#print "randx=",x

			y=(2*random.random()-1.)*maxL
			a=(2*random.random()-1.)*maxL
			#self.positions[i]=(x,y,a)
			self.positions[i]=(x,y,a)


	def moveGrain(self,i,dx,dy,da):
		
		angle=self.positions[i][2]
		cx=self.positions[i][0]
		cy=self.positions[i][1]
	
		#if dx>self.L:
		#	print "WOAH"
	        # update positions
		cx=cx+dx
		cy=cy+dy
		angle=angle+da

	
	        # go periodic 
		cx=(cx+self.L/2.)%self.L-self.L/2.
		cy=(cy+self.L/2.)%self.L-self.L/2.

		"""
                
		if (cx>self.L/2.):
			cx=cx-self.L
		elif (cx<-self.L/2.):
			cx=cx+self.L
		if (cy>self.L/.2):
			cy=cy-self.L
		elif (cy<-self.L/2.):
			cy=cy+self.L
                
		"""

                #print self.positions

	        #update positions
		#self.positions[i][0]=cx
		#self.positions[i][1]=cy
		#self.positions[i][2]=angle

                self.positions[i]=(cx,cy,angle)

	def testOverlapByIndex(self,i,j):
		pi=int(self.shapeTypes[i])
		pj=int(self.shapeTypes[j])

		xi=self.positions[i][0]
		yi=self.positions[i][1]
		ai=self.positions[i][2]

		xj=self.positions[j][0]
		yj=self.positions[j][1]
		aj=self.positions[j][2]


		anyOverlap=0
		self.movePolygon(pi,xi,yi,ai)
		self.movePolygonJ(pj,xj,yj,aj)

		doesOverlap=self.polygons[pi].overlaps(self.polygonsJ[pj])
		if doesOverlap==True:
			anyOverlap=anyOverlap+1


		# periodic
		if (xi-xj) > (self.L/2.):
			xj=xj+self.L
		if (xi-xj) < (-self.L/2.):
			xj=xj-self.L
		if (yi-yj) > (self.L/2.):
			yj=yj+self.L
		if (yi-yj) < (-self.L/2.):
			yj=yj-self.L

		# create a polygon modeled on the original shape
		
		self.movePolygon(pi,xi,yi,ai)
		self.movePolygonJ(pj,xj,yj,aj)

		doesOverlap=self.polygons[pi].overlaps(self.polygonsJ[pj])

		if doesOverlap==True:
			anyOverlap=anyOverlap+1
		
		if anyOverlap>0:
			return True
		else:
			return False

	def findR2Cutoff(self,polygonInput):
		b=polygonInput.boundingBox()
		r2max=0.
		for i in range(len(b)):
			for j in range(len(b)):
				if i!=j:
					r2=(b[i]-b[j])**2
					if r2>r2max:
						r2max=r2
		return r2max

	def findR2CutoffMax(self):
		r2CutoffMax=0.
		for pi in range(len(self.polygons)):
			#thisPoly=self.polygons[pi]
			self.movePolygon(pi,0.,0.,0.)
			thisPoly=self.polygons[pi]
			thiscut=self.findR2Cutoff(thisPoly)
			if thiscut>r2CutoffMax:
				r2CutoffMax=thiscut
		#print "find found:", r2CutoffMax
		return r2CutoffMax
			

        def findCenter(self,polygonInput):

            b=polygonInput.boundingBox()
            xmin=b[0]
            xmax=b[1]
            ymin=b[2]
            ymax=b[3]
            
            cx=(xmax-xmin)/2.+xmin
            cy=(ymax-ymin)/2.+ymin

            return (cx,cy)

        def drawBoundingBox(self,polygonInput):
            b=polygonInput.boundingBox()
            xmin=b[0]
            xmax=b[1]
            ymin=b[2]
            ymax=b[3]

            curve(pos=[(xmin,ymax,0.),(xmax,ymax,0.),(xmax,ymin,0.),(xmin,ymin,0.),(xmin,ymax,0.)])

            curve(pos=[(xmin,ymax,0.),(xmax,ymax,0.)],color=color.blue)

            print "xmin=",xmin
            print "xmax=",xmax

            cx=(xmax-xmin)/2.+xmin
            cy=(ymax-ymin)/2.+ymin
            sphere(pos=(cx,cy,0.),radius=.2)
            
	
	def setupViaData(self,input_data):

		#FILE=open(filename,'r')

	
		inPickle=input_data
			

		self.initialContours=inPickle[1]
		self.shapeTypes=inPickle[2]
		self.L=inPickle[3]
		self.N=inPickle[4]
		self.pressure=inPickle[5]
		self.volume=inPickle[6]
		self.replicaLabel=inPickle[7]
		algorithmDetails=inPickle[8]
		self.t=inPickle[9]
		#self.initialSize=inPickle[10]
		finalSize=inPickle[10]


		self.polygons=[]
		self.polygonsJ=[]
		#self.polygonsDisplay=[]
		self.polygonPositions=[]
		self.polygonPositionsJ=[]
		

		# convert position values to tuples:
		for p in range(len(inPickle[0])):
			thisPos=inPickle[0][p]
			#print thisPos
			#print "thisPos=",thisPos
			thisPos=tuple(thisPos)
	
			self.positions.append(thisPos)
		

		#self.positions=inPickle[0]

		#print "self.positions=",self.positions

		# fill out the polygons
		r2Max=0.

		for i in range(len(self.initialContours)):

			thisPolygon=Polygon.Polygon(self.initialContours[i])
			thisPolygonJ=Polygon.Polygon(self.initialContours[i])
			      
			#print "test of contour=",thisPolygon

			"""
			thisPolygon=Polygon.Polygon(self.polygonContoursInitial[i])
			self.boundingBoxesInitial.append(thisPolygon.boundingBox())
			thisPolygonJ=Polygon.Polygon(self.polygonContoursInitial[i])
			#thisPolygonDisplay=Polygon.Polygon(self.polygonContoursInitial[i])
			"""
			
		     
			self.polygons.append(thisPolygon)
			self.polygonsJ.append(thisPolygonJ)
			#self.polygonsDisplay.append(thisPolygonDisplay)
			


			self.polygonPositions.append((0.,0.,0.))
			self.polygonPositionsJ.append((0.,0.,0.))
			#self.polygonPositionsDisplay.append((0.,0.,0.))
		

		
		# scale the polygons

		#initialSize=self.polygons[0].boundingBox()[0]
		initialSize=abs(self.polygons[0].boundingBox()[1]-self.polygons[0].boundingBox()[0])
		scaleRatio=float(finalSize)/float(initialSize)
		
		#box_i=self.boundingBoxesInitial[0]
		#box_f=self.boundingBoxesFinal[0]
		
		#sideInitial=box_i[1]-box_i[0]
		#sideFinal=box_f[1]-box_f[0]
		#scaleRatio=sideFinal/float(sideInitial)

		for pi in range(len(self.polygons)):
			(cx,cy)=self.polygons[pi].center()
			self.polygons[pi].scale(scaleRatio,scaleRatio,cx,cy)
			(cx,cy)=self.polygonsJ[pi].center()
			self.polygonsJ[pi].scale(scaleRatio,scaleRatio,cx,cy)
                

		self.r2cutoff=self.findR2CutoffMax()
		
		self.updateNeighborList(self.addedR2Buffer)

		#print "r2cutoff=",self.r2cutoff

		#get the density
		totalParticleArea=0.
		for i in range(self.N):
			thisShape=int(self.shapeTypes[i])
			totalParticleArea=totalParticleArea+self.polygons[thisShape].area()
		totalParticleArea=self.getTotalParticleArea()

		self.density=totalParticleArea/self.L**2
		self.volume=self.N/self.density
		
		# set the initial translate step stats
		
		self.translate_step=algorithmDetails[0]
		self.rotate_step=algorithmDetails[1]
		self.scale_step=algorithmDetails[2]


		"""
		if len(algorithmDetails)>0:
			print "yes"
			
		else:

			self.translate_step=sqrt(self.polygons[0].area())/16.
			self.rotate_step=math.pi/16.
			self.scale_step=.001*sqrt(self.polygons[0].area())
		"""

		
        def setupTest(self):

           
            self.N=12
            self.L=15.

 
            for i in range(self.N):    
                self.shapeTypes.append(0)
                self.positions.append((2*i-self.L/2.,0.,0.))
                
            
            self.N=len(self.positions)

            cx=0.
            cy=0.
            
            R=.5

            self.translateStep=R/3.
            self.rotateStep=math.pi/2.

            theta=math.pi/3.

            mid=(R/2.)*math.tan(theta/2.)/3.

            thisTri=[(cx-R/2.,cy-R*math.sin(theta)/2.),(cx+R/2.,cy-R*math.sin(theta)/2.),(cx,cy+R*math.sin(theta)/2.)]

	    #thisTri=[(cx-R/2.,cy+R),(cx,cy+R*1.5),(cx+R/2.,cy+R),(cx+R/2.,cy-R),(cx,cy-R*1.5),(cx-R/2.,cy-R)]

	    #thisTri=[(cx-R/2.,cy+R/2.),(cx,cy+R*1.5),(cx+R/2.,cy+R),(cx+R/2.,cy-R),(cx,cy-R*1.5),(cx-R/2.,cy-R)]

	    thisTri=[(cx-R/2.,cy+R/2.),(cx+R/2.,cy+R/2.),(cx+R/2.,cy-R/2.),(cx-R/2.,cy-R/2.)]

	    #thisTri=[(cx-R/2.,cy-R*math.sin(theta)/2.),(cx+R/2.,cy-R*math.sin(theta)/2.),(cx,cy+R*math.sin(theta)/2.),(cx-R,cy+R*math.sin(theta)/2.+R)]
	    
	    polygonContours=[]

            polygonContours.append(thisTri)
            thisPolygon=Polygon.Polygon(polygonContours[0])
            thisPolygonJ=Polygon.Polygon(polygonContours[0])
            #thisPolygonDisplay=Polygon.Polygon(polygonContours[0])
            self.polygons.append(thisPolygon)
            self.polygonsJ.append(thisPolygonJ)
            #self.polygonsDisplay.append(thisPolygonDisplay)
            self.polygonPositions.append((0.,0.,0.))
            self.polygonPositionsJ.append((0.,0.,0.))
            #self.polygonPositionsDisplay.append((0.,0.,0.))

            center=thisPolygon.center()
            print "center=",center, "mid=",mid
            
	    

	    self.translate_step=sqrt(self.polygons[0].area())/16.
	    self.rotate_step=math.pi/16.
	    self.scale_step=.001*sqrt(self.polygons[0].area())
 
        def moveGrainSweep(self):

		#self.assignCells()
		#self.findNeighborsR2()
		#print self.r2neighborlist

		for i in range(self.N):

			
			randi=random.randint(0,self.N)
			

			# decide whether this is a translation or rotation move
			moveType=random.randint(0,2)
			#moveiype=1

			#print "moveType=",moveType

			
			dx=0.
			dy=0.
			da=0.

			if moveType==0:

				dx=2*random.random()*self.translateStep-self.translateStep
				dy=2*random.random()*self.translateStep-self.translateStep

				da=0.
			if moveType==1:
				dx=0.
				dy=0.
				da=2*random.random()*self.rotateStep-self.rotateStep


			self.moveGrain(randi,dx,dy,da)

			#print self.findAnyOverlapsWithNeighbors(randi)
			#time.sleep(1)


			
			if (self.findAnyOverlapsWithNeighbors(randi)==True):
				self.moveGrain(randi,-dx,-dy,-da)
				print "movedBack"

			else:
				if moveType==0:
					self.translateAttemptSuccess=self.translateAttemptSuccess+1
				if moveType==1:
					self.rotateAttemptSuccess=self.rotateAttemptSuccess+1

			if moveType==0:
				self.translateAttemptCounter=self.translateAttemptCounter+1
			if moveType==1:
				self.rotateAttemptCounter=self.rotateAttemptCounter+1
	

        def runRand(self,runTime):
		for i in range(0,runTime):

			for i in range(0,self.N):
				#print i
				randi=random.randint(0,self.N-1)
				dx=(2*random.random()-1.)*.01
				dy=(2*random.random()-1.)*.01
				da=(2*random.random()-1.)*.01
				self.moveGrain(randi,dx,dy,da)
			self.t=self.t+1
		
		
		

	def run(self,runTimeInSweeps):

		
		self.dVscale=0.
		self.dVabsolute=0.

	        self.thisTime=0. 
		self.thisSweep=0.
		#print self.density, self.positions[0]
		#print "self.thisSweep=",self.thisSweep
		#print "thisSweep=",self.thisSweep,"runTimeInSweeps=",runTimeInSweeps
		while (self.thisSweep<runTimeInSweeps):

	            #MOVE THIS DOWN!!!!
                    self.t=self.t+1
		    #print "self.t=",self.t
		    #self.sweeps=float(self.t/float(self.N))
		    #print "sweeps=",self.sweeps
		    #self.sweeps=float(self.t/float(self.N))
		    #self.thisTime=self.thisTime+1
		    #self.thisSweep=self.thisTime/float(self.N) 

		    #update the neighborList
		    if (self.t%self.r2neighborUpdateCount==0):
			    self.updateNeighborList(self.addedR2Buffer)
			    #print "neighbor example:",self.r2neighborlist[0]

		    moveChooser=random.random()

		    if moveChooser<self.p_translate:
			moveType=0
		    if moveChooser>=self.p_translate and moveChooser < self.p_rotate:
			moveType=1
		    if moveChooser>=self.p_rotate:
			moveType=2

		    #moveType=random.randint(0,3)

		    #randi=random.randint(0,self.N-1)

		    dx=0.
		    dy=0.
		    da=0.

		    if moveType==0:

			randi=random.randint(0,self.N-1)
			dx=(2*random.random()-1.)*self.translate_step
			dy=(2*random.random()-1.)*self.translate_step
			da=0.
			self.moveGrain(randi,dx,dy,da)

			self.translate_attempt=self.translate_attempt+1

		    if moveType==1:
			randi=random.randint(0,self.N-1)
			dx=0.
			dy=0.
			da=(2*random.random()-1.)*self.rotate_step
			self.moveGrain(randi,dx,dy,da)

			self.rotate_attempt=self.rotate_attempt+1


		    if moveType==2:

			boltzReject=0

			self.dVscale=(2*random.random()-1.)*self.scale_step+1. #hovers around 1
			self.dVabsolute=self.dVscale*self.L**2-self.L**2
			#print dVabsolute

			# reject any move that isn't allowed vi boltzmann statistics
			if self.dVabsolute>0.:
			    self.boltzFactor=math.exp(-self.pressure*(self.dVabsolute))
			    self.dVAcceptProb=min([1,self.boltzFactor])
			    #print self.pressure, dVabsolute, dVAcceptProb

			    if random.random()<self.dVAcceptProb:
				self.scalePolygons(1/self.dVscale) # 1 over the volume expansion
				self.scale_attempt=self.scale_attempt+1
			    else:
				boltzReject=1

			else:
			    #print self.pressure, dVabsolute
			    self.scalePolygons(1/self.dVscale) # 1 over the volume expansion
			    self.r2cutoff=self.findR2CutoffMax()
			    self.scale_attempt=self.scale_attempt+1



		    if moveType==2:

			if boltzReject==0:

			    anyOverlap=0
			    for i in range(0,self.N): # note: going over *all* particles
				for j in range(0,self.N):
				    if j!=i:
					doesOverlap=self.testOverlapByIndex(i,j)
					if (doesOverlap==True):
					    anyOverlap=anyOverlap+1
			    if anyOverlap>0:
				#print "overlap"
				if moveType==0 or moveType==1:
				    self.moveGrain(randi,-dx,-dy,-da)
				if moveType==2:
				    self.scalePolygons(self.dVscale)
				    # update the r2cutoff
				    self.r2cutoff=self.findR2CutoffMax()
				    #print "r2cutoff B =",self.r2cutoff
			    else:
				#update the density and the volume
				#self.density=self.getTotalParticleArea()/self.L**2
				#self.volume=self.N/self.density # assume that the area of all original particles is such that N*(area of all original particles summed)=1
				self.scale_accept=self.scale_accept+1


		    if moveType==1 or moveType==0:
			# update sweeps
			self.thisSweep=self.thisSweep+float(1./self.N)
			self.sweeps=self.sweeps+float(1./self.N)

			anyOverlap=0
			i=randi # note: only testing randi
			for j in range(0,self.N):
			    if j!=i:
				doesOverlap=self.testOverlapByIndex(i,j)
				if (doesOverlap==True):
				    anyOverlap=anyOverlap+1
			if anyOverlap>0:
			    #print "overlap"
			    if moveType==0 or moveType==1:
				self.moveGrain(randi,-dx,-dy,-da)
			    if moveType==2:
				self.scalePolygons(1/dV)
				# update the r2cutoff
				self.r2cutoff=self.findR2CutoffMax()
				#print "r2cutoff A =",self.r2cutoff 
			else:
			    #update the sweeps:
			    #self.thisSweep=self.thisSweep+float(1./self.N)
                            #self.sweeps=self.sweeps+float(1./self.N)
			    if moveType==0:
                                #print "dx,dy=",dx,dy
				dr=math.sqrt(dx**2+dy**2)
				
				self.ave_velocity=self.ave_velocity+dr
				self.total_displacement=self.total_displacement+dr
				self.translate_accept=self.translate_accept+1
			    if moveType==1:
				self.rotate_accept=self.rotate_accept+1

		#update density and volume -- NEW, 8/12/11
		self.density=self.getTotalParticleArea()/self.L**2
                self.volume=self.N/self.density		

		#test of the particle coords:
		#print "inside runner:"
		#print " ... after loop:"
		#print self.density, self.positions[0]

		#update statistics
		if (self.t%self.stat_collect==0):


			#print "translate_step/L=",self.translate_step/self.L

			#print "t=",self.t,"updating"
			if self.rotate_attempt>0:
			    self.rotate_accept=self.rotate_accept/float(self.rotate_attempt)

			    # modify the step
			    epsilon=-self.accept_goal+self.rotate_accept
			    self.rotate_step=self.rotate_step+self.K_c*epsilon*self.rotate_step
			    self.rotate_step=self.rotate_step%(3.14159*2.)

			    self.rotate_attempt=0

			#print "self.translate_step=",self.translate_step
			#print "self.denisty=",self.density

			if self.translate_attempt>0:

			
			    if (self.translate_accept>0):
				    #self.ave_velocity=self.ave_velocity/float(self.translate_accept)
				    self.ave_velocity=self.total_displacement/float(self.translate_accept)
				    #print "ave_velocity=",self.ave_velocity
				    #print "total displ=",self.total_displacement
				    #print "trans accept=",self.translate_accept
				    #print "r2cutff=",self.r2cutoff
				    #self.ave_velocity=self.total_displacement/float(self.translate_accept)
				    self.total_displacement=0.
				    if (self.ave_velocity>0):
					    timeTakenToMoveOneParticleSpan=math.sqrt(self.r2cutoff)/self.ave_velocity
					    self.r2cutoffUpdateCount=int(math.ceil(.5*timeTakenToMoveOneParticleSpan))
					    #print "new count=",self.r2cutoffUpdateCount

			    #self.r2cutoffUpdateCount=int(.5*math.ceil(math.sqrt(self.r2cutoff)/self.ave_velocity))
			    #print "updateCount=",self.r2cutoffUpdateCount
				    #self.ave_velocity=0.

			    self.translate_accept=self.translate_accept/float(self.translate_attempt)

			    #print "trans accept ratio=",self.translate_accept
			    # modify the step
			    epsilon=-self.accept_goal+self.translate_accept
			    self.translate_step=self.rotate_step+self.K_c*epsilon*self.translate_step
			    self.translate_step=self.translate_step%self.L # to make sure we don't go too big

			    self.translate_attempt=0

			if self.scale_attempt>0:
			    self.scale_accept=self.scale_accept/float(self.scale_attempt)

			    # modify the step
			    epsilon=-self.accept_goal+self.scale_accept
			    self.scale_step=self.scale_step+self.K_c*epsilon*self.scale_step


			    self.scale_attempt=0



	def saveProgressFile(self,fileLabel):

		progressVals=str(self.t)+" "+str(self.N)+" "+str(self.density)+" "+str(self.pressure)+" " + str(self.volume)

		

		progressName="N_"+str(self.N)+"_"+fileLabel+"_progress.out"

		progressLine=progressVals+"\n"

		FILE=open(progressName,'a')
		FILE.write(progressLine)
		FILE.close()

		return progressVals

	def generatePicklePackage(self):
	
		self.movePolygon(0,0.,0.,0.)

		# calculate the density right now:
		#totalParticleArea=self.getTotalParticleArea()
		#self.density=totalParticleArea/self.L**2
		#self.volume=self.N/self.density
		
		finalSize=abs(self.polygons[0].boundingBox()[1]-self.polygons[0].boundingBox()[0])
		
		algorithmDetails=[self.translate_step,self.rotate_step,self.scale_step]

		picklePackage=[self.positions,self.initialContours,self.shapeTypes,self.L,self.N,self.pressure,self.volume,self.replicaLabel,algorithmDetails,self.t,finalSize,self.density]

		return picklePackage

 
	def saveLatestSnapshot(self,fileLabel):
		snapFileName="N_"+str(self.N)+"_"+fileLabel+"_snap.out"
		#snapFileName="N_"+str(self.N)+"_replicaLabel_"+str(self.replicaLabel)+"_snap.out"
		FILE=open(snapFileName,'w')
		
		#print "fileLabel=",fileLabel
		#print "len(self.polygonContoursInitial)=",len(self.polygonContoursInitial)
		# get the contours
		#polygonContours=[]

		# get the final contours
		currentContours=[]
		for i in range(len(self.polygons)):
			thisPolygon=self.polygons[i]
			self.movePolygon(i,0.,0.,0.)
			currentContours.append(thisPolygon[0])
	
		algorithmDetails=[self.translate_step,self.rotate_step,self.scale_step]

		picklePackage=[self.positions,currentContours,self.shapeTypes,self.L,self.N,self.pressure,self.volume,self.replicaLabel,algorithmDetails,self.t]

		#systemProperties=[self.L,self.N,self.pressure,self.volume,self.replicaLabel]
		#algorithmProperties=[self.simTime,self.translateStep,self.rotateStep,self.volumeStep,self.K_c_move,self.moveStepAssessInterval,self.targetMoveAcceptanceRatio,self.r2cutoff]
		#picklePackage=[grainProperties,systemProperties]
		pickle.dump(picklePackage,FILE)
		FILE.close()
		
		#print self.translate_step, self.rotate_step,self.scale_step

		

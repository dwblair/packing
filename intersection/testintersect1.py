import numpy as np
#from scipy import *
from pylab import *

def clip(subjectPolygon, clipPolygon):
   def inside(p):
      return(cp2[0]-cp1[0])*(p[1]-cp1[1]) > (cp2[1]-cp1[1])*(p[0]-cp1[0])
 
   def computeIntersection():
      dc = np.array([ cp1[0] - cp2[0], cp1[1] - cp2[1] ])
      dp = np.array([ s[0] - e[0], s[1] - e[1] ])
      n1 = cp1[0] * cp2[1] - cp1[1] * cp2[0]
      n2 = s[0] * e[1] - s[1] * e[0] 
      n3 = 1.0 / (dc[0] * dp[1] - dc[1] * dp[0])
      return np.array([(n1*dp[0] - n2*dc[0]) * n3, (n1*dp[1] - n2*dc[1]) * n3])
 
   outputList = np.copy(subjectPolygon)
   cp1 = clipPolygon[-1]
 
   for clipVertex in clipPolygon:
      cp2 = clipVertex
      inputList = np.copy(outputList)
      outputList = []
      print inputList
      s = inputList[-1]
 
      for subjectVertex in inputList:
         e = subjectVertex
         if inside(e):
            if not inside(s):
               outputlist=np.append(outputList,computeIntersection())
            outputList=np.append(outputList,e)
         elif inside(s):
             outputlist=np.append(outputList,computeIntersection())
         s = e
      cp1 = cp2
   return(outputList)

#def drawPoly(poly):


   
poly1=[[0,0],[10,0],[10,10],[0,10],[0,0]]
poly2=[[5,0],[10,0],[10,20],[5,20],[0,0]]
poly3=clip(poly1,poly2)

a=np.array(poly3)
print a
a=np.array(poly1)
b=np.array(poly2)
c=np.array(poly3)
plot(a[:,0],a[:,1],label='a')
plot(b[:,0],b[:,1],label='b')
plot(c[:,0],c[:,1],label='c')
legend()
show()


import scipy
import sys
import matplotlib.pyplot as plt

#inputfilelist=['uv_1.0.txt','uv_1.1.txt','uv_1.5.txt','uv_1.9.txt','uv_2.3.txt','uv_2.7.txt','uv_3.1.txt','uv_3.5.txt']
#aspectRatios=[1.,1.1,1.5,1.9,2.3,2.7,3.1,3.5]

inputfilelist=['uv_3.5.txt']
aspectRatios=[3.5]

#inputfile=sys.argv[1]
for i in range(0,len(aspectRatios)):
#for inputfile in inputfilelist:
    inputfile=inputfilelist[i]
    f=open(inputfile,'r')
    lines=f.readlines()
    f.close
#lines=lines[len(lines)-51:len(lines)+1]
    u=[]
    v=[]
    for line in lines:
        these=line.split(',')
        u.append(float(these[0]))
        v.append(float(these[1]))
    uMin=min(u)
    uflip=[]
    vflip=[]
    #print uMin
    #u=u-uMin
    for j in range(0,len(u)):
        uflip.append(-1.*float(u[j]))
        vflip.append(-1.*float(v[j]))
    #for j in range(0,len(u)):
    #    u[j]=u[j]-uMin
    plt.plot(u,v,'o',label='A='+str(aspectRatios[i]))
    #uflip=-u
    #vflip=-v
    #plt.plot(vflip,vflip)
plt.xlim(-4.,4.)
plt.ylim(-4.,4.)
plt.legend(loc='lower left')
plt.show()

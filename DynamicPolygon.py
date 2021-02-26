import numpy as np
import math as m 
import matplotlib.pyplot as plt

#A particle front approach to edge detection for 
#composite yarn cross-section outline
#taken from μCT data

def mag(V):
    return np.sqrt(V.dot(V))
def diff(V0,V1):
	dV=V0-V1
	return np.sqrt(dV.dot(dV))
class FixedNode:
    def __init__(self,PositionV,Value):
       self.PositionV=PositionV
       self.Value=Value

class Node:
    def __init__(self,Index,RefNode,PositionV,Speed0,Value,Iterations):
        self.Index=Index
        self.PositionV=PositionV
        DPosVector=PositionV-RefNode.PositionV 
        UDPosVector=DPosVector/mag(DPosVector)        
        self.SpeedV=UDPosVector*Speed0
        self.SpeedV0=UDPosVector*Speed0
        self.Value=Value
        self.DValue0=abs(Value-RefNode.Value)
        self.Control=False
        self.Activations=[]
        self.Trace=np.empty([Iterations,2])
        self.FactorThreshold=10.0
        self.Factors=[]
    def UpdatePosition(self,Dt):
        NewPos=self.PositionV+self.SpeedV*Dt
        self.PositionV=NewPos#np.array([int(NewPos[0]),int(NewPos[1])])
    def UpdateSpeed(self,RefNode,Value,ReverseSpeedFactor):
        DValue=abs(Value-RefNode.Value)
        self.Value=Value
        Factor=abs(DValue-self.DValue0)
        #self.DValue0=DValue
        self.Factors.append(Factor)
        #self.FactorThreshold=Factor*1.5
        self.FactorThreshold=(sum(self.Factors)/len(self.Factors))
        #print(Factor)
        aa=25 # Activations allowed
        check=Factor-self.FactorThreshold
        if check>10.0:
          self.Activations.append(Factor)
          if len(self.Activations)>aa: 
             self.Control=True	
             self.SpeedV=-ReverseSpeedFactor*self.SpeedV
          else:
             self.SpeedV=(1/check)*self.SpeedV   
        #else:     
         # self.SpeedV=self.SpeedV0

if __name__=='__main__':
    	
  path='D:\\IJPythonReconstructionOfTexComp\\'
  imgName='X498'   #image name  	
  GSA=np.genfromtxt(path+imgName+'.txt')
  Dt=1.0 #Time step
  iters=300 # Solution iterations
  Speed0=0.5 #Initial particle speed 
  Num=100 # Number of particles
  w0=12 #Center Pixel window size
  w1=8# Particle Pixel window size
  ReverseSpeedFactor=0.1# Corrective backstep speed factor (if needed)
  CVec=np.array([304,376])#[237,168])#Centre point (fixed) taken from click on ImageJ
  CValue=np.average(GSA[CVec[0]-w0:CVec[0]+w0,CVec[1]-w0:CVec[1]+w0]) # Compute base average intensity
  C=FixedNode(CVec,CValue)
  R=3 #radius of polygon initialisation circle
  Dang=(2*m.pi)/(Num-1)# Angle step
  NodeList=[]
  #Initilise output plot
  img = plt.imread(path+imgName+'.png')
  
  fig, ax = plt.subplots()
  
  ax.imshow(img)
  ax.scatter(CVec[1],CVec[0])
  ax.annotate('Centre',(CVec[1],CVec[0]))
  #Initialise polygon nodes
  for i in range(Num):
    pos=np.array([m.cos(Dang*i)*R+CVec[0],m.sin(Dang*i)*R*4.0+CVec[1]])
    value=np.average(GSA[int(pos[0]-w1):int(pos[0]+w1),int(pos[1]-w1):int(pos[1]+w1)])
    N=Node(i,C,pos,Speed0,value,iters)
    NodeList.append(N)
  #Run simulation
  for N in NodeList:
    upto=iters-1
    for i in range(iters):
       N.UpdatePosition(Dt)
       pos=N.PositionV
       N.Trace[i]=pos
       if N.Control:
          upto=i
          break
       X=int(pos[0])
       Y=int(pos[1])
       NewPVal=np.average(GSA[X-w1:X+w1,Y-w1:Y+w1])
       N.UpdateSpeed(C,NewPVal,ReverseSpeedFactor)

    #Add nodes to plot
    #ax.scatter(N.Trace[0,1],N.Trace[0,0])
    ax.scatter(N.Trace[upto,1],N.Trace[upto,0])
    #for p in [0,upto]:
    ax.annotate(str(N.Index),(N.Trace[upto,1],N.Trace[upto,0]))
    print(N.Index,len(N.Activations))
  plt.show()
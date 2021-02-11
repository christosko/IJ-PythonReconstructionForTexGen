import numpy as np
import matplotlib.pyplot as plt
#import marplotlib.image as img
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
        self.Trace=np.empty([Iterations,2])
    def UpdatePosition(self,Dt):
        NewPos=self.PositionV+self.SpeedV*Dt
        self.PositionV=NewPos
    def UpdateSpeed(self,RefNode,Value,FactorThreshold,ReverseSpeedFactor):
        DValue=abs(Value-RefNode.Value)
        Factor=DValue/self.DValue0
        print(Factor)
        if Factor>FactorThreshold:
          self.Control=True	
          self.SpeedV=-ReverseSpeedFactor*self.SpeedV
        else:     
          self.SpeedV=self.SpeedV*0.95	

if __name__=='__main__':
  # For full polygon: 
  #1)Determine number of points and radius
  #2)Compute position vectors     	
  path='D:\\IJPythonReconstructionOfTexComp\\'
  imgName='X430'        	
  GSA=np.genfromtxt(path+imgName+'.txt')
  Dt=1.0
  iters=100
  Speed0=5
  #Initialise
  
  r=3
  FactorThreshold=6.0
  ReverseSpeedFactor=0.5
  CPVec=np.array([237,168])
  PPVec0=np.array([224,190])
  CValue=np.average(GSA[CPVec[0]-r:CPVec[0]+r,CPVec[1]-r:CPVec[1]+r])
  PValue0=np.average(GSA[PPVec0[0]-r:PPVec0[0]+r,PPVec0[1]-r:PPVec0[1]+r])
  C=FixedNode(CPVec,CValue)
  P=Node(0,C,PPVec0,Speed0,PValue0)  
  for i in range(iters):
     P.UpdatePosition(Dt)
     pos=P.PositionV
     P.Trace[i]=pos
     if P.Control:
        upto=i+1
        break
     X=int(pos[0])
     Y=int(pos[1])
     NewPVal=np.average(GSA[X-r:X+r,Y-r:Y+r])
     P.UpdateSpeed(C,NewPVal,FactorThreshold,ReverseSpeedFactor)
     
  img = plt.imread(path+imgName+'.png')
  
  fig, ax = plt.subplots()
  
  ax.imshow(img)
  ax.scatter(CPVec[1],CPVec[0])
  ax.annotate('Centre',(CPVec[1],CPVec[0]))
  ax.scatter(P.Trace[:upto,1],P.Trace[:upto,0])
  for i in range(upto):
    ax.annotate(str(i),(P.Trace[i,1],P.Trace[i,0]))

  plt.show()
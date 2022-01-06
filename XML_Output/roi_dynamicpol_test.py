from ij import IJ
from ij.gui import PointRoi, PolygonRoi, Roi
import math
from jarray import zeros, array
from Jama import Matrix
# Get current ImagePlus
image = IJ.getImage()
# Get current ROI
roi = image.getRoi()
slic=image.getSlice()
### Define : 
insert='X_0_0_'
####
#A particle front approach to edge detection for 
#composite yarn cross-section outline
#taken from μCT data

def mag(V):
    return ((V.times(V.transpose())).get(0,0))**0.5
def diff(V0,V1):
	dV=V0.minus(V1)
	return ((dV.times(dV.transpose())).get(0,0))**0.5
def FL(N):
    return math.floor(N)
def CL(N):
    return math.ceil(N)
def CO(N):
    return math.cos(N)
def SI(N):
    return math.sin(N)
                	
def Ellipse(R,cx,cy,fx,fy,n,Dang,i):
    return Matrix([[CO(Dang*i)*R*fx+cx, (FL(SI(Dang*i))+CL(SI(Dang*i)))*(((FL(SI(Dang*i))+CL(SI(Dang*i)))*SI(Dang*i))**n)*R*fy+cy]])            
def average(i,j,w,img):
  sum=0
  i=int(i)
  j=int(j)
  w=int(w)
  for x in range(-w,w):
    for y in range(-w,w):
       Val=img.getPixel(i+x,j+y)
       sum+=Val[0]       
  return sum/(2*w+1)**2

class FixedNode:
    def __init__(self,PositionV,Value):
       self.PositionV=PositionV
       self.Value=Value       
class Node:
    def __init__(self,Index,RefNode,PositionV,Speed0,Value,Iterations):
        self.Index=Index
        self.PositionV=PositionV
        DPosVector=PositionV.minus(RefNode.PositionV) 
        UDPosVector=DPosVector.timesEquals(1/mag(DPosVector))        
        self.SpeedV=UDPosVector.timesEquals(Speed0)
        self.SpeedV0=UDPosVector.timesEquals(Speed0)
        self.Value=Value
        self.DValue0=abs(Value-RefNode.Value)
        self.Control=False
        self.Activations=[]
        #self.Trace=Matrix(Iterations,2)
        self.FactorThreshold=15.0
        self.Factors=[]
    def UpdatePosition(self,Dt):
        NewPos=self.PositionV.plus(self.SpeedV.timesEquals(Dt))
        self.PositionV=NewPos
    def UpdateSpeed(self,RefNode,Value,ReverseSpeedFactor,aa,Dth):
        DValue=abs(Value-RefNode.Value)
        self.Value=Value
        Factor=abs(DValue-self.DValue0)
        #self.DValue0=DValue
        self.Factors.append(Factor)
        self.FactorThreshold=(sum(self.Factors)/len(self.Factors))
        check=Factor-self.FactorThreshold
        if check>Dth:
          self.Activations.append(Factor)
          if len(self.Activations)>aa: 
             self.Control=True	
             self.SpeedV=self.SpeedV.timesEquals(-ReverseSpeedFactor)
          else:
             self.SpeedV=self.SpeedV.timesEquals(1/check)  
Dt=1.0 #Time step
iters=1000 # Solution iterations
Speed0=0.5 #Initial particle speed 
Num=65 # Number of particles
ReverseSpeedFactor=0.1# Corrective backstep speed factor (if needed) 
Dang=(2*math.pi)/(Num-1)# Angle step
aa=100 # Activations allowed
Dth=8.0 # Difference threshold
NodeList=[]

#w0 Center Pixel window size
#w1  Particle Pixel window size
#R  Radius of initial polygon
#fx Ellipse horizontal stretch 
#fy  Ellipse vertical stretch 
#n  Power ellipse

mode='yarn'
if mode in ['Yarn','yarn']:
  w0,w1,R,fx,fy,n=4,4,4,3.0,0.8,0.2
elif mode in ['Void', 'void']:
  w0,w1,R,fx,fy,n=1,1,1,1.0,1.0,1.0


if roi:
    # Get ROI points
    polygon = roi.getPolygon()
    n_points = polygon.npoints
    x = polygon.xpoints
    y = polygon.ypoints
    val=image.getPixel(x[0],y[0])
    CVec=Matrix([[x[0],y[0]]])#[237,168])#Centre point (fixed) taken from user point
    CValue=average(x[0],y[0],w0,image) # Compute base average intensity
    C=FixedNode(CVec,CValue)             

    for i in range(Num):
       cx=x[0]
       cy=y[0]

       pos=Ellipse(R,cx,cy,fx,fy,n,Dang,i)
       #pos=Matrix([[math.cos(Dang*i)*R*fx+cx, math.sin(Dang*i)*R*fy+cy]])
       value=average(pos.get(0,0),pos.get(0,1),w1,image)
       N=Node(i,C,pos,Speed0,value,iters)
       NodeList.append(N)
   
    for N in NodeList:
       upto=iters-1
       for i in range(iters):
          N.UpdatePosition(Dt)
          pos=N.PositionV
          #N.Trace.setMatrix(pos)
          if N.Control:
            upto=i
            break
          X=int(pos.get(0,0))
          Y=int(pos.get(0,1))
          NewPVal=average(X,Y,w1,image)
          N.UpdateSpeed(C,NewPVal,ReverseSpeedFactor,aa,Dth) 
    PRoi=PointRoi()
    for N in NodeList:
       Pos=N.PositionV.getArray()[0]
       PRoi.addPoint(Pos[0],Pos[1])
    image.setRoi(PRoi,0)
    image.updateImage()
    
Save=0     
dirPath="D:\\IJPythonReconstructionOfTexComp\\Data2\\"
if Save:
   fileName=insert+str(slic)+".dat"
   f = open(dirPath+fileName, 'w')
   for N in NodeList:
      Pos=N.PositionV.getArray()[0]  
      #print str(Pos[0]) + " " + str(Pos[1])
      f.write(str(Pos[0]) + " " + str(Pos[1]) + "\n")
   f.close()
   print "Saved: " + fileName

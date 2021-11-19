#2D Structure tensor orientation
from ij import IJ
from ij.gui import PointRoi, PolygonRoi, Roi
import math as m
from jarray import zeros, array
from Jama import Matrix, EigenvalueDecomposition
import imp
import os
import sys

sys.path.append(os.getcwd())
os.chdir(os.getcwd())
#import StructureTensorDecomp as st 

#imp.reload(StructureTensorDecomp)

def value(i,j,img):
  return img.getPixel(i,j)[0] 

def averagedIdx2(i,j,h,w,img):
  sum=0.0
  h=int(h)
  den=12*h
  i=int(i)
  j=int(j)
  w=int(w)
  for x in range(-w,w):
    for y in range(-w,w):
       Ix_2h=value(i+x-2*h,j+y,img)
       Ix_h=value(i+x-h,j+y,img)
       Ixh=value(i+x+h,j+y,img)
       Ix2h=value(i+x+2*h,j+y,img)
       dIdx=(Ix_2h-8*Ix_h+8*Ixh-Ix2h)/den      
       sum+=dIdx**2      
  return sum/(2*w+1)**2
  
def averagedIdy2(i,j,h,w,img):
  sum=0.0
  h=int(h)
  den=12*h
  i=int(i)
  j=int(j)
  w=int(w)
  for x in range(-w,w):
    for y in range(-w,w):
       Iy_2h=value(i+x,j+y-2*h,img)
       Iy_h=value(i+x,j+y-h,img)
       Iyh=value(i+x,j+y+h,img)
       Iy2h=value(i+x,j+y+2*h,img)
       dIdy=(Iy_2h-8*Iy_h+8*Iyh-Iy2h)/den      
       sum+=dIdy**2      
  return sum/(2*w+1)**2

def averagedIdxdIdy(i,j,h,w,img):
  sum=0.0
  h=int(h)
  den=12*h
  i=int(i)
  j=int(j)
  w=int(w)
  for x in range(-w,w):
    for y in range(-w,w):
       Ix_2h=value(i+x-2*h,j+y,img)
       Ix_h=value(i+x-h,j+y,img)
       Ixh=value(i+x+h,j+y,img)
       Ix2h=value(i+x+2*h,j+y,img)
       dIdx=(Ix_2h-8*Ix_h+8*Ixh-Ix2h)/den   
           
       Iy_2h=value(i+x,j+y-2*h,img)
       Iy_h=value(i+x,j+y-h,img)
       Iyh=value(i+x,j+y+h,img)
       Iy2h=value(i+x,j+y+2*h,img)
       dIdy=(Iy_2h-8*Iy_h+8*Iyh-Iy2h)/den      
       sum+=dIdx*dIdy     
  return sum/(2*w+1)**2

def StructureTensor(i,j,h,w,img):
  h=int(h)
  w=int(w)
  den=float(12.0*h)
  i=int(i)
  j=int(j)
  
  h=float(h)

  ST11=averagedIdx2(i,j,h,w,img)#dIdx**2
  ST22=averagedIdy2(i,j,h,w,img)#dIdy**2
  #ST12=dIdx*dIdy #Needs to be replaced with bivariate 2nd order derivative
  ST12=averagedIdxdIdy(i,j,h,w,img)
  ST=Matrix([[ST11, ST12],[ST12,ST22]])
  EDST=EigenvalueDecomposition(ST)
  Eval=EDST.getRealEigenvalues()
  EV=EDST.getV()

  return Eval, EV
# Get current ImagePlus
image = IJ.getImage()


roi = image.getRoi()

steps=100
d=8
w=12
h=2

polygon = roi.getPolygon()
n_points = polygon.npoints
x = polygon.xpoints
y = polygon.ypoints

pos=Matrix([[x[n_points-1],y[n_points-1]]])
PRoi=PointRoi()
for i in range(steps):
   eigval,eigvec=StructureTensor(pos.get(0,0),pos.get(0,1),h,w,image)
   scaling=5
   Min=20000.000
   minind=0
   for e in range(len(eigval)):
      if eigval[e]<Min:
         Min=eigval[e]
         minind=e
   Max=0.0000
   for e in range(len(eigval)):
      if eigval[e]>Max:
         Max=eigval[e]
      if Max>0:
         ratio=1-Min/Max
      else:
         ratio=0.0
     #print(ratio)  
   scaling=scaling*ratio
   evec0=eigvec.getMatrix(minind,minind,0,1)
   pos=pos.plus(evec0.times(2.0))
   x=pos.get(0,0)
   y=pos.get(0,1)
   PRoi.addPoint(x, y)

image.setRoi(PRoi,0)
image.updateImage()   


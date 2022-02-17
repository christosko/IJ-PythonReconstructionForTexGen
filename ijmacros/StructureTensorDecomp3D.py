#3D Structure tensor orientation
from ij import IJ, ImagePlus 
from ij.gui import PointRoi, PolygonRoi, Roi
from ij.gui import Line, Overlay
import math as m
from jarray import zeros, array
from Jama import Matrix, EigenvalueDecomposition
# Get current ImagePlus
imp=IJ.getImage()
image = imp.getImageStack()
#iproc=process.ImageProcessor()
# Get current ROI
#roi = image.getRoi()
#slic=image.getSlice()
def average(i,j,k,w,img):
  sum=0.0
  i=int(i)
  j=int(j)
  k=int(k)
  w=int(w)
  for x in range(-w,w):
    for y in range(-w,w):
      for z in range(-w,w):
         Val=img.getVoxel(i+x,j+y,k+z)
         sum+=Val[0]       
  return sum/(2*w+1)**3
  
def value(i,j,k,img):
  return img.getVoxel(i,j,k)
   
def averagedIdx(i,j,k,h,w,img):
  sum=0.0
  h=int(h)
  den=12*h
  i=int(i)
  j=int(j)
  k=int(k)
  w=int(w)
  for x in range(-w,w):
    for y in range(-w,w):
      for z in range(-w,w):
       Ix_2h=value(i+x-2*h,j+y,k+z,img)
       Ix_h=value(i+x-h,j+y,k+z,img)
       Ixh=value(i+x+h,j+y,k+z,img)
       Ix2h=value(i+x+2*h,j+y,k+z,img)
       dIdx=(Ix_2h-8*Ix_h+8*Ixh-Ix2h)/den      
       sum+=dIdx      
  return sum/(2*w+1)**3
  
def averagedIdy(i,j,k,h,w,img):
  sum=0.0
  h=int(h)
  den=12*h
  i=int(i)
  j=int(j)
  w=int(w)
  k=int(k)
  for x in range(-w,w):
    for y in range(-w,w):
      for z in range(-w,w):
        Iy_2h=value(i+x,j+y-2*h,k+z,img)
        Iy_h=value(i+x,j+y-h,k+z,img)
        Iyh=value(i+x,j+y+h,k+z,img)
        Iy2h=value(i+x,j+y+2*h,k+z,img)
        dIdy=(Iy_2h-8*Iy_h+8*Iyh-Iy2h)/den      
        sum+=dIdy      
  return sum/(2*w+1)**3

def averagedIdz(i,j,k,h,w,img):
  sum=0.0
  h=int(h)
  den=12*h
  i=int(i)
  j=int(j)
  w=int(w)
  k=int(k)
  for x in range(-w,w):
    for y in range(-w,w):
      for z in range(-w,w):
        Iz_2h=value(i+x,j+y,k+z-2*h,img)
        Iz_h=value(i+x,j+y,k+z-h,img)
        Izh=value(i+x,j+y,k+z+h,img)
        Iz2h=value(i+x,j+y,k+z+2*h,img)
        dIdz=(Iz_2h-8*Iz_h+8*Izh-Iz2h)/den      
        sum+=dIdz      
  return sum/(2*w+1)**3  

def averagedIdx2(i,j,k,h,w,img):
  sum=0.0
  h=int(h)
  den=12*h
  i=int(i)
  j=int(j)
  k=int(k)
  w=int(w)
  for x in range(-w,w):
    for y in range(-w,w):
      for z in range(-w,w):
        Ix_2h=value(i+x-2*h,j+y,k+z,img)
        Ix_h=value(i+x-h,j+y,k+z,img)
        Ixh=value(i+x+h,j+y,k+z,img)
        Ix2h=value(i+x+2*h,j+y,k+z,img)
        dIdx=(Ix_2h-8*Ix_h+8*Ixh-Ix2h)/den      
        sum+=dIdx**2      
  return sum/(2*w+1)**3
  
def averagedIdy2(i,j,k,h,w,img):
  sum=0.0
  h=int(h)
  den=12*h
  i=int(i)
  j=int(j)
  k=int(k)
  w=int(w)
  for x in range(-w,w):
    for y in range(-w,w):
      for z in range(-w,w):
        Iy_2h=value(i+x,j+y-2*h,k+z,img)
        Iy_h=value(i+x,j+y-h,k+z,img)
        Iyh=value(i+x,j+y+h,k+z,img)
        Iy2h=value(i+x,j+y+2*h,k+z,img)
        dIdy=(Iy_2h-8*Iy_h+8*Iyh-Iy2h)/den      
        sum+=dIdy**2      
  return sum/(2*w+1)**3

def averagedIdz2(i,j,k,h,w,img):
  sum=0.0
  h=int(h)
  den=12*h
  i=int(i)
  j=int(j)
  w=int(w)
  k=int(k)
  for x in range(-w,w):
    for y in range(-w,w):
      for z in range(-w,w):
        Iz_2h=value(i+x,j+y,k+z-2*h,img)
        Iz_h=value(i+x,j+y,k+z-h,img)
        Izh=value(i+x,j+y,k+z+h,img)
        Iz2h=value(i+x,j+y,k+z+2*h,img)
        dIdz=(Iz_2h-8*Iz_h+8*Izh-Iz2h)/den      
        sum+=dIdz**2      
  return sum/(2*w+1)**3  

def averagedIdxdIdy(i,j,k,h,w,img):
  sum=0.0
  h=int(h)
  den=12*h
  i=int(i)
  j=int(j)
  k=int(k)
  w=int(w)
  for x in range(-w,w):
    for y in range(-w,w):
      for z in range(-w,w):
        Ix_2h=value(i+x-2*h,j+y,k+z,img)
        Ix_h=value(i+x-h,j+y,k+z,img)
        Ixh=value(i+x+h,j+y,k+z,img)
        Ix2h=value(i+x+2*h,j+y,k+z,img)
        dIdx=(Ix_2h-8*Ix_h+8*Ixh-Ix2h)/den   
            
        Iy_2h=value(i+x,j+y-2*h,k+z,img)
        Iy_h=value(i+x,j+y-h,k+z,img)
        Iyh=value(i+x,j+y+h,k+z,img)
        Iy2h=value(i+x,j+y+2*h,k+z,img)
        dIdy=(Iy_2h-8*Iy_h+8*Iyh-Iy2h)/den      
        sum+=dIdx*dIdy     
  return sum/(2*w+1)**3

def averagedIdxdIdz(i,j,k,h,w,img):
  sum=0.0
  h=int(h)
  den=12*h
  i=int(i)
  j=int(j)
  k=int(k)
  w=int(w)
  for x in range(-w,w):
    for y in range(-w,w):
      for z in range(-w,w):
        Ix_2h=value(i+x-2*h,j+y,k+z,img)
        Ix_h=value(i+x-h,j+y,k+z,img)
        Ixh=value(i+x+h,j+y,k+z,img)
        Ix2h=value(i+x+2*h,j+y,k+z,img)
        dIdx=(Ix_2h-8*Ix_h+8*Ixh-Ix2h)/den   
            
        Iz_2h=value(i+x,j+y,k+z-2*h,img)
        Iz_h=value(i+x,j+y,k+z-h,img)
        Izh=value(i+x,j+y,k+z+h,img)
        Iz2h=value(i+x,j+y,k+z+2*h,img)
        dIdz=(Iz_2h-8*Iz_h+8*Izh-Iz2h)/den      
        sum+=dIdx*dIdz     
  return sum/(2*w+1)**3

def averagedIdydIdz(i,j,k,h,w,img):
  sum=0.0
  h=int(h)
  den=12*h
  i=int(i)
  j=int(j)
  k=int(k)
  w=int(w)
  for x in range(-w,w):
    for y in range(-w,w):
      for z in range(-w,w):
        Iy_2h=value(i+x,j+y-2*h,k+z,img)
        Iy_h=value(i+x,j+y-h,k+z,img)
        Iyh=value(i+x,j+y+h,k+z,img)
        Iy2h=value(i+x,j+y+2*h,k+z,img)
        dIdy=(Iy_2h-8*Iy_h+8*Iyh-Iy2h)/den 
            
        Iz_2h=value(i+x,j+y,k+z-2*h,img)
        Iz_h=value(i+x,j+y,k+z-h,img)
        Izh=value(i+x,j+y,k+z+h,img)
        Iz2h=value(i+x,j+y,k+z+2*h,img)
        dIdz=(Iz_2h-8*Iz_h+8*Izh-Iz2h)/den      
        sum+=dIdy*dIdz     
  return sum/(2*w+1)**3

def StructureTensor(i,j,k,h,w,img):
  h=int(h)
  w=int(w)
  den=float(12.0*h)
  i=int(i)
  j=int(j)
  k=int(k)
  
  h=float(h)

  ST11=averagedIdx2(i,j,k,h,w,img)#dIdx**2
  ST22=averagedIdy2(i,j,k,h,w,img)#dIdy**2
  ST33=averagedIdz2(i,j,k,h,w,img)
  #ST12=dIdx*dIdy #Needs to be replaced with bivariate 2nd order derivative
  ST12=averagedIdxdIdy(i,j,k,h,w,img)
  ST13=averagedIdxdIdz(i,j,k,h,w,img)
  ST23=averagedIdydIdz(i,j,k,h,w,img)  
  ST=Matrix([[ST11, ST12, ST13],[ST12, ST22, ST23],[ST13, ST23, ST33]])
  EDST=EigenvalueDecomposition(ST)
  Eval=EDST.getRealEigenvalues()
  EV=EDST.getV()

  return Eval, EV

d=48
w=4
h=2

Roi=imp.getRoi()

fbounds=Roi.getFloatBounds()

grid_size=5
xrng=range(int(fbounds.x),int(fbounds.x+fbounds.width),grid_size)
yrng=range(int(fbounds.y),int(fbounds.y+fbounds.height),grid_size)
S=imp.getSlice()

zrng=range(int(S),int(S+1))

#PRoi=PointRoi()
Over=Overlay()

for i in xrng:
  for j in yrng:
    for k in zrng:
      scaling=15
      eigval,eigvec=StructureTensor(i,j,k,h,w,image)  
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
      #print(Min,Max)  
      scaling=scaling*ratio
      x=i
      y=j
      z=int(k)
      if ratio>0.7:
         #print(x,y,z)               
         # print(ratio,eigvec.get(minind,0),eigvec.get(minind,1),eigvec.get(minind,2))
         
         #PRoi.addPoint(x,y)  
         x1=float(x+eigvec.get(minind,0)*scaling)
         y1=float(y+eigvec.get(minind,2)*scaling)
         #z1=int(z+eigvec.get(minind,1)*(scaling/2))
         L=Line.create(x,y,x1,y1)
         L.setWidth(2)
         L.setPosition(zrng[0])
         Over.add(L)
         #PRoi.addPoint(x1,y1)
         #x2=float(x+eigvec.get(minind,0)*scaling)
         #y2=float(y+eigvec.get(minind,1)*scaling)
         #z2=int(z+eigvec.get(minind,1)*scaling)
         #PRoi.addPoint(x2,y2) 
      else:
         pass
         #PRoi.addPoint(x,y)  
#imp.setRoi(PRoi,0)
#imp.updateImage()
imp.setOverlay(Over)
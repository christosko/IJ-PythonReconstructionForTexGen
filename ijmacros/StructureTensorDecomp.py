#2D Structure tensor orientation
from ij import IJ
from ij.gui import PointRoi, PolygonRoi, Roi
from ij.gui import Line, Overlay
import math as m
from jarray import zeros, array
from Jama import Matrix, EigenvalueDecomposition
# Get current ImagePlus
image = IJ.getImage()
# Get current ROI
#roi = image.getRoi()
#slic=image.getSlice()
def average(i,j,w,img):
  sum=0.0
  i=int(i)
  j=int(j)
  w=int(w)
  for x in range(-w,w):
    for y in range(-w,w):
       Val=img.getPixel(i+x,j+y)
       sum+=Val[0]       
  return sum/(2*w+1)**2
  
def value(i,j,img):
  return img.getPixel(i,j)[0] 
   
def averagedIdx(i,j,h,w,img):
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
       sum+=dIdx      
  return sum/(2*w+1)**2
  
def averagedIdy(i,j,h,w,img):
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
       sum+=dIdy      
  return sum/(2*w+1)**2

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

d=6
w=18
h=2
#PRoi=PointRoi()
Over=Overlay()
Roi=image.getRoi()

fbounds=Roi.getFloatBounds()

grid_size=5

xrng=range(int(fbounds.x),int(fbounds.x+fbounds.width),grid_size)
yrng=range(int(fbounds.y),int(fbounds.y+fbounds.height),grid_size)

for i in xrng:
  for j in yrng:
     eigval,eigvec=StructureTensor(i,j,h,w,image)
     scaling=15
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
     x=i
     y=j
     if ratio>0.999:               
       x1=float(x+eigvec.get(minind,0)*scaling)
       y1=float(y-eigvec.get(minind,1)*scaling)
       #z1=int(z+eigvec.get(minind,1)*(scaling/2))
       L=Line.create(x,y,x1,y1)
       L.setWidth(2)
       #L.setPosition(zrng[0])
       Over.add(L)     
       #PRoi.addPoint(i,j)  
       #PRoi.addPoint(float(i+eigvec.get(minind,0)*(scaling/2)),float(j+eigvec.get(minind,1)*(scaling/2)))
       #PRoi.addPoint(float(i+eigvec.get(minind,0)*scaling),float(j+eigvec.get(minind,1)*scaling)) 
     else:
       pass
       #PRoi.addPoint(i,j)  
#image.setRoi(PRoi,0)
image.setOverlay(Over)
image.updateImage()
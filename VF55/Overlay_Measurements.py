#3D import all measured yarn sections and overlay them 
from ij import IJ, ImagePlus
from ij.gui import PointRoi, PolygonRoi, Roi, Overlay
import math as m
from jarray import zeros, array
from Jama import Matrix, EigenvalueDecomposition
import os 
import sys
from os import listdir

# Get current ImagePlus
imp=IJ.getImage()
image = imp.getImageStack()
# Get measured polygons
cwd=os.getcwd()
dataloc=cwd+'//VF55//Data7'
dats=[d for d in listdir(dataloc) if '.dat' in d and 'X' in d] # Change for slice direction
##
MyOverlay=Overlay()

for i,d in enumerate(dats):
  info=(d.replace('.dat','')).split('_')
  Direction=file[0]
  YarnIndex=int(file[1])
  Index=int(file[2])
  Slice=int(file[3])
  Polygon=PointRoi() 
  fil=open(d,'r')
  lines=fil.readlines()  
  sumx=0
  sumy=0
  N=0
  for l in lines:
    coo=l.split(' ')
    x=int(coo[0])
    y=int(coo[1])
    Polygon.addPoint(x,y)
    sumx+=x
    sumy+=y
    N+=1
  ax=sumx//N
  ay=sumy//N  
  Polygon.setPosition(Slice)
  MyOverlay.add(Polygon)
  Label=TextRoi(ax,ay,str(YarnIndex),8)
  Label.setPosition(Slice)
  MyOverlay.add(Label)
image.setOverlay(MyOverlay)
 
  

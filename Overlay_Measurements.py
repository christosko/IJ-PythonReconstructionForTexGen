#3D import all measured yarn sections and overlay them 
from ij import IJ, ImagePlus
from ij.gui import PointRoi, PolygonRoi, Roi, Overlay, TextRoi
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
dataloc='D:\\IJPythonReconstructionOfTexComp\\VF55\\Data8'
#Axis:
ax='Y'
dats=[d for d in listdir(dataloc) if '.dat' in d and ax in d] # Change for slice direction
##
MyOverlay=Overlay()

for i,d in enumerate(dats):
  info=(d.replace('.dat','')).split('_')
  Direction=info[0]
  YarnIndex=int(info[1])
  Index=int(info[2])
  Slice=int(info[3])
  Polygon=PointRoi() 
  print d
  fil=open(dataloc+'\\'+d,'r')
  lines=fil.readlines()  
  sumx=0
  sumy=0
  N=0
  for l in lines:
    coo=l.split(' ')
    x=int(coo[0])
    y=int(coo[2])
    Polygon.addPoint(x,y)
    sumx+=x
    sumy+=y
    N+=1
  ax=sumx//N
  ay=sumy//N  
  Polygon.setPosition(Slice)
  MyOverlay.add(Polygon)
  Label=TextRoi(ax,ay,str(YarnIndex))
  Label.setPosition(Slice)
  MyOverlay.add(Label)
imp.setOverlay(MyOverlay)
 
  

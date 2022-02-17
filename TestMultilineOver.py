#Tests overlay of lines
from ij import IJ, ImagePlus
from ij.gui import PointRoi, PolygonRoi, Roi
from ij.gui import Line, Overlay
import math as m
import random as ra

imp=IJ.getImage()
image = imp.getImageStack()

Over=Overlay()
x=range(200,300,5)
y=range(200,300,5)

#PRoi.addPoint(x,y)
for px in x:
  for py in y:
      extx=ra.randrange(-15, 15)
      exty=ra.randrange(-15, 15)
      px1=px+extx
      py1=py+exty
      L=Line.create(px,py,px1,py1)
      L.setWidth(2)
      L.setPosition(276)
      Over.add(L)
imp.setOverlay(Over)         
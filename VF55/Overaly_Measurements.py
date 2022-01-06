#3D Structure tensor orientation
from ij import IJ, ImagePlus
from ij.gui import PointRoi, PolygonRoi, Roi
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
dats=[d.replace('.dat','') for d in listdir(dataloc) if '.dat' in d]


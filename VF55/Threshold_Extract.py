from ij import IJ, ImagePlus
from ij.gui import PointRoi, PolygonRoi, Roi
import math as m
from jarray import zeros, array
from Jama import Matrix, EigenvalueDecomposition
import os 
# Get current ImagePlus
imp=IJ.getImage()
image = imp.getImageStack()
#
#Set voxel size in mm:
svox=0.0105

#
workdir=os.getcwd()
os.chdir(workdir)
dim=imp.getDimensions()
thRngVoid=range(0,147)
thRngMatrix=range(147,210)
thRngFibre=range(210,255)

VoidData=open('VoidXYZ.dat','w')
MatrixData=open('MatrixXYZ.dat','w')
FibreData=open('FibreXYZ.dat','w')

for i in range(dim[0]):
  for j in range(dim[1]):
    for k in range(dim[3]):
      gsval=image.getVoxel(i,j,k)
      if gsval in thRngVoid:
        VoidData.write(str(i*svox*0.5)+' '+str(j*svox*0.5)+' '+str(k*svox*0.5)+'\n')
      elif gsval in thRngMatrix:
        MatrixData.write(str(i*svox*0.5)+' '+str(j*svox*0.5)+' '+str(k*svox*0.5)+'\n')
      elif gsval in thRngFibre:
        FibreData.write(str(i*svox*0.5)+' '+str(j*svox*0.5)+' '+str(k*svox*0.5)+'\n') 
      else:
        pass
VoidData.close()
MatrixData.close()
FibreData.close()             
      
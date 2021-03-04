from ij import IJ
from ij.gui import PolygonRoi, Roi
import math
#import os 
#import sys------------------------------------------------------------------------------------------------
#index=input('Give yarn index:')
#import csv
#yarn_indx=0
#yarn_type=''
# Get current ImagePlus
image = IJ.getImage()

# Get current ROI
roi = image.getRoi()
slic=image.getSlice()
if roi:
    # Get ROI points
    polygon = roi.getPolygon()
    n_points = polygon.npoints
    x = polygon.xpoints
    y = polygon.ypoints

    dirPath='D:\\IJPythonReconstructionOfTexComp\\Data\\'
   # os.chdir(dirPath)
    fileName='Z_4_1_'+str(slic)+'.dat'

    f = open(dirPath+fileName, 'w')
    try:
        #writer = csv.writer(f)
        #writer.writerow(('x', 'y'))
        for i in range(n_points):
            f.write(str(x[i])+'  '+str(y[i])+'\n')
    finally:
        f.close()
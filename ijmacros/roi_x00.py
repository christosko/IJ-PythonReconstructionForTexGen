from ij import IJ
from ij.gui import PolygonRoi, Roi
import math
import os
# Get current ImagePlus
image = IJ.getImage()
# Get current ROI
roi = image.getRoi()
slic=image.getSlice()
### Define : 
insert='Y_1_0_'
###

if roi:
    # Get ROI points
    polygon = roi.getPolygon()
    n_points = polygon.npoints
    x = polygon.xpoints
    y = polygon.ypoints

    dirPath='D:\\IJPythonReconstructionOfTexComp\\Data2\\'
    fileName=insert+str(slic)+'.dat'
    if os.path.isdir(dirPath):
      f = open(dirPath+fileName, 'w')
    else: 
      try:
        os.mkdir(dirPath)
        f = open(dirPath+fileName, 'w')
      except OSError: 
        print 'Not authorised to create:'+dirPath
    
    try:
        for i in range(n_points):
            f.write(str(x[i])+'  '+str(y[i])+'\n')
    finally:
        f.close()
    print 'Saved: '+fileName
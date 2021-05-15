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
insert='X_33_0_' + str(slic) # + '_-1'
###

if roi:
    # Get ROI points
    polygon = roi.getPolygon()
    n_points = polygon.npoints
    x = polygon.xpoints
    y = polygon.ypoints

    dirPath='D:\\IJPythonReconstructionOfTexComp\\VF55\\Data3\\'
    fileName=insert + '.dat'
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
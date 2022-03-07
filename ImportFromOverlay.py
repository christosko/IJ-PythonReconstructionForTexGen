from ij import IJ, ImagePlus
from ij.gui import PolygonRoi, Roi
import math
import os
# Get current ImagePlus
image = IJ.getImage()
# Get current ROI
roi = image.getRoi()
#roi=roi.getInterpolatedPolygon(5,True)
slic=image.getSlice()
### Define : 

MyOver=image.getOverlay()
I=0
roi=MyOver.get(I)   
while roi:
    # Get ROI points
    #polygon = roi.getPolygon()
    points=roi.getContainedPoints()
    print(len(points))
    pos=roi.getPosition()
    print('slice:',pos)
    #print(points[0].x,points[0].y)
    print(roi.getTypeAsString())
    textroi=MyOver.get(I+1)

    insert='X_'+(textroi.getText()).replace('\n','')+'_0_' + str(pos)# + '_-1'
    dirPath='D:\\IJPythonReconstructionOfTexComp\\VF64\\Data3\\'
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
        for p in points:
            f.write(str(p.x)+'  '+str(p.y)+'\n')
    finally:
        f.close()
    print 'Saved: '+fileName
    I+=2
    roi=MyOver.get(I)
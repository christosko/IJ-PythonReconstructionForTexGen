from ij import IJ
from ij.gui import PolygonRoi, Roi
import math
from jarray import zeros, array
from Jama import Matrix
# Get current ImagePlus
image = IJ.getImage()
# Get current ROI
roi = image.getRoi()
slic=image.getSlice()
myx=array([142,243,82],'i')
M=Matrix([[3,4]])

print M.times(M.transpose()).get(0,0)
if roi:
    # Get ROI points
    polygon = roi.getPolygon()
    n_points = polygon.npoints
    x = polygon.xpoints
    y = polygon.ypoints
    val=image.getPixel(x[0],y[0])
    print val
    print x, y   

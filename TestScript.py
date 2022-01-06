import numpy as np
import math as m
import os
import sys
from os import listdir
import matplotlib.pyplot as plt
cwd=os.getcwd()
sys.path.append(cwd)
#Requires a minimum build of texgen and a path pointing to Python libraries
TexGenPath='C:\\Python27\\Lib\\site-packages\\TexGen'
sys.path.append(TexGenPath)

from TexGen.Core import*

tosort=XYZVector()
tosort.push_back(XYZ(2.3,4.5,3.2))
tosort.push_back(XYZ(4.3,3.1,3.2))
tosort.push_back(XYZ(33.1,2.7,3.2))
tosort.push_back(XYZ(1.3,10.5,3.2))
tosort.push_back(XYZ(21.3,9.2,3.2))

print(tosort)
vecsort=sorted(tosort,key=[0].y)
print(vecsort)
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

p=XYZ(23.03,433.04,5.4)
print(p*float(1.0/50.0))

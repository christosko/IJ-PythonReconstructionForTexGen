import os,sys
import math

import platform
from os import listdir
import matplotlib.pyplot as plt
import numpy as np
from math import*
#from numpy import *

cwd=os.getcwd()

DatName='VectorField_Y_Fine.txt'

fig1, ax1=plt.subplots()

ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_title("Coherency 0.2 - 0.4")


Table=np.genfromtxt(DatName)
X=Table[:,0]
Y=Table[:,1]
Z=Table[:,2]
Coh=Table[:,6]

#plt.hist(Coh, bins='auto')
#plt.show()
  

Sub=[(i,val) for i,val in enumerate(Coh) if val>0.2 and val<0.4 ]



SubX=[X[i[0]] for i in Sub if Z[i[0]]==276]
SubY=[-Y[i[0]] for i in Sub if Z[i[0]]==276]

ax1.scatter(SubX,SubY)
plt.show()
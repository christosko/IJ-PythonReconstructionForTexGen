import numpy as np
import math as m
import os
import sys
from os import listdir
cwd=os.getcwd()
sys.path.append(cwd)
#Requires a minimum build of texgen and a path pointing to Python libraris
TexGenPath='C:\\Python27\\Lib\\site-packages\\TexGen'
sys.path.append(TexGenPath)

from TexGen.Core import*


def Extend(Textile):
    Yarns=Textile.GetYarns()
    Weft=[]
    Warp=[]
    Binder=[]
    Xunit=XYZ(1.0,0.0,0.0)
    Yunit=XYZ(0.0,1.0,0.0)

    for i,Y in enumerate(Yarns):
       
       M0P=Y.GetNode(0).GetPosition()
       M1P=Y.GetNode(1).GetPosition()
       V=M1P-M0P
       M2=Y.GetNode(2)
       if M2:
          print 'Binder : '+str(i)
          Binder.append(Y)
       elif GetLength(CrossProduct(Xunit,V))<GetLength(CrossProduct(Yunit,V)):
          print 'Warp : '+str(i)
          Warp.append(Y)
       elif GetLength(CrossProduct(Xunit,V))>GetLength(CrossProduct(Yunit,V)):  
          print 'Warp : '+str(i)
          Weft.append(Y)
           

    return NewTextile
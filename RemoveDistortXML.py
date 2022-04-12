import numpy as np
import math as m
import os
import sys
from os import listdir

import xml.etree.ElementTree as ET
from xml.dom import minidom

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D 
from sklearn.cluster import KMeans

import imp
import ControlMesh as cm

imp.reload(cm)
cwd=os.getcwd()
sys.path.append(cwd)
#Requires a minimum build of texgen and a path pointing to Python libraries
TexGenPath='C:\\Python27\\Lib\\site-packages\\TexGen'
sys.path.append(TexGenPath)

from TexGen.Core import*

def GetDomain(Path,TG3Name):
   Tree=ET.parse(Path+'\\'+TG3Name)
   Root=Tree.getroot()
   domain=[]
   for child in Root:
     Textile=child
     break
   for child in Textile:
     if child.tag=='Domain':
        Domain=child
        for child in Domain:
           if child.tag=='Plane':
              domain.append(float(child.attrib['d']))
   CP2=XYZ(domain[0],domain[2],domain[4])
   CP1=XYZ(-domain[1],-domain[3],-domain[5])
   CDomain=CDomainPlanes(CP1,CP2)           
   return CDomain 
   
def ImportFromTG3(ModelLocation,ModelName,TexName):
   ReadFromXML(ModelLocation+'\\'+ModelName+'.tg3')
   Textile=GetTextile(TexName)
   return Textile


if __name__=='__main__':
  TestPath='D:\\Project_Repository\\Damage_Analysis\\S11'
  TGName='OrthRec8thSheared.tg3'
  TGTexName='ShearedDomain'
  Textile=ImportFromTG3(TestPath,TGName,TGTexName)
  OldDomain=GetDomain(TestPath,TGName)
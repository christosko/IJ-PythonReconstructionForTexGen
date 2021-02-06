
import xml.etree.ElementTree as ET
from xml.dom import minidom
import numpy as np
#from ij import IJ
from ij.gui import PolygonRoi, Roi
import math
import os
#import sys
#import csv
#Define XML formating function, sorting functions, data structures 
def PrintXML(Root,Path,Name):
  rough_string =ET.tostring(Root, 'utf-8')
  reparsed = minidom.parseString(rough_string)
  Out=reparsed.toprettyxml(indent="  ")
  F=open(Path+xmlName,'w')
  F.writelines(Out)
  F.close()  
  return 0
def SortFiles(FileNames):
 switch=True
 while switch:
    sw=False
    for i in range(len(FileNames)-1):
       f=int(FileNames[i].split('_')[2])
       fn=int(FileNames[i+1].split('_')[2])
       if int(f)>int(fn):
          temp=FileNames[i]
          FileNames[i]=FileNames[i+1]
          FileNames[i+1]=temp
          sw=True
    switch=sw
 return FileNames 

class Yarn:

   def __init__(self,Index):
      self.Index=Index
      self.Section=None
      self.left=None
      self.right=None

   def InsertSection(self,Index,Section):

      if self.Index:

         if Index<self.Index:
           if self.left is None:
             self.left=Yarn(Index)
           else:
             self.left.InsertSection(Index,Section)

         elif Index>self.Index:

           if self.right is None:
             self.right=Node(Index)
           else:
             self.right.InsertSection(value)

      else:
         self.Index=Index

   def find(self,inval):
      if inval<self.value:
         if self.left is None:
            return str(inval)+'Not Found'
         return self.left.find(inval)
      elif inval>self.value:
         if self.right is None:
            return str(inval)+'Not Found'
         return self.right.find(inval)
      else:
         print(str(self.value)+'is found')

   def PrintTree(self):
      if self.left:
        self.left.PrintTree()
      print(self.value)
      if self.right:
        self.right.PrintTree()

class Section:
  def __init__(self,t,Polygon)
    self.t=t
    self.Polygon=Polygon
    self.left=None
    self.right=None

  def InsertPolygon(self,t,Polygon):
    return  
#Initialise an XML File with the basic root and children

resolution=0.0105
gsize=np.array([901 611 532])
InName='TestXML_2.xml'
Root=ET.Element('XML',{'version':'1.0'})
TGModel=ET.SubElement(Root,'TexGenModel',{'version':'3.12.0'})
TextileProps={

 'name':"G",
 'GeometryScale':"mm", 
 'type':"CTextile", 
 'NeedsBuilding':"0"
  
}
Textile=ET.SubElement(TGModel,'Textile',TextileProps)

Path='D:\\Polygon_Data_55VF\\'
os.chdir(Path)
#Find window size in pixels: 
#file=open('window_size.txt','r')
#WinSize=file.read()
#file.close()
WinSize=np.genfromtxt('window_size.txt')

#Find pixel size: 
file=open('pixel_size.txt','r')
ImgRes=file.read() #mm
file.close()
ImgRes=float(ImgRes)
#Define domain
DS=WinSize*ImgRes

NormalsStr=[
"1.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00",
"-1.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00",
"0.000000000000e+00, 1.000000000000e+00, 0.000000000000e+00",
"0.000000000000e+00, -1.000000000000e+00, 0.000000000000e+00",
"0.000000000000e+00, 0.000000000000e+00, 1.000000000000e+00",
"0.000000000000e+00, 0.000000000000e+00, -1.000000000000e+00"
]

D=[
str(0.0),
str(-DS[0]),
str(0.0),
str(-DS[1]),
str(-DS[2]),
str(0.0)
]

Domain=SubElement(Textile,'Domain',{'type':'CDomainPlanes'})

Planes=[ET.Element('Plane', {'Normal':NormalsStr[i],'d':D[i]})]
Domain.Extend(Planes)

YarnDict={
'YarnLinearDensityValue':"6.700000000000e-05" 
'YarnLinearDensityUnits':"kg/m" 
'FibreDensityValue':"1.720000000000e+00" 
'FibreDensityUnits':"g/cm^3" 
'NumSlaveNodes':"50" 
'NumSectionPoints':"100" 
'NeedsBuilding':"7"
}

WarpYarns={}
WeftYarns={}
BinderYarns={}  

FileList=[f for f in listdir(Path)]
FileList=[(f.replace('.dat','')).split('_') for f in listdir(Path)]



X=[[F for F in FileList if F[0]=='X'] for i in ]
Y=[F for F in FileList if F[0]=='Y']
BZ=[F for f in FileList if F[0]=='BZ']
BX=[F for f in FileList if F[0]=='BX'] 



for FN in FileNames:
  FN=FN.replace('.dat','')
  FN=FN.split('_')
  Type=FN[0]
  Index=FN[1]
  Slice=int(FN[2])

  Types.append(Type)
  Indices.append(Index)
  Slices.append(Slice)

LastYarnI=int(max(Indices))
for y in Range(LastYarnI+1):
  YFiles=[FileNames[i] for i in range(FileNames) if Indices[i]==y]
  SortedFiles=SortFiles(YFiles)
  YType=Types[y]
  for f in SortedFiles:

  for i in range(len(SortedFiles)):
    YarnDict[Slice[]]=np.genfromtxt(f)

for i in  
 YarnSorted=[]


#you need first and last slice. we need to decide wether we need common slices across yarns


#InName='TestXML.xml'
#OutName='OutXML.xml'
#Tree=ET.parse(Path+xmlName)
#Root=Tree.getroot()

#for child in Root:
#    Textile=child
#    break
for i in range(len(Slices)):

  first_slice=min([Slices[j] for j in range(len(Slices)) if Indices[j]==Indices[i]])
  last_slice=max([Slices[j] for j in range(len(Slices)) if Indices[j]==Indices[i]])
  t=(Slices[i]-first_slice)/(last_slice-first_slice)

  Yarn_Found=False
  for child in Textile:
    if child.tag=='Yarn':
      Yarns_Found+=1
      index=child.attrib['index']
      if index==Indicies[i]:
        yarn=child
        Yarn_Found=True
        break

  points=np.genfromtxt(FileNames[i])
  points=points*ImgRes
  x=[points[k][0] for k in range(len(points))]
  y=[points[k][1] for k in range(len(points))]
  x0=sum(x)/len(x)
  y0=sum(y)/len(y)

  if not Yarn_Found:

    yarn=ET.SubElement(Textile,'Yarn',{'index':Indices[i]})
    interpolation=ET.SubElement(Yarn,'Interpolation',{'Periodic':1,
      'ForceInPlaneTagent':0,'type':'CInterpolationBezier'})
    YarnSection==ET.SubElement(Yarn,'YarnSection',{'type':'CYarnSectionInterpPosition',
     'ConstMesh':1, 'Ramped':1, 'Polar':0})
    ind=0
  else:
    ind=1


  PositionSection=ET.SubElement(YarnSection,'PositionSection',{'t':str(t)})
  Section=ET.SubElement(PositionSection,'Section', {'type':'CSectionPolygon'})
  SectionMesh=ET.SubElement(Section,'SectionMesh',{'type':'CSectionMeshRectangular',
   'NumLayer':'2', 'TriangleCorners':'1'})

  PolygonPoints=[ET.Element('PolygonPoint',{'value':str(x[j])+', '+str(y[j])}) for j in range(len(points))]
  Section.extend(PolygonPoints)
  FibreDistribution=ET.SubElement(yarn,'FibreDistribution',{'type':'CFibreDistributionConst'})
  if Types[i]=='X':
    MasterNode0=ET.SubElement(yarn,'MasterNode', {'index':str(ind),
       'Position': str(slices[i]*ImgRes)+', '+str(x0)+', '+str(y0),
       'Tangent':'0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00',
       'Up':'0.000000000000e+00, 0.000000000000e+00, 1.000000000000e+00',
       'Angle':'0.0' })
  elif Types[i]=='Y':
    MasterNode0=ET.SubElement(yarn,'MasterNode', {'index':str(ind),
       'Position': str(x0)+', '+str(slices[i]*ImgRes)+', '+str(y0),
       'Tangent':'0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00',
       'Up':'0.000000000000e+00, 0.000000000000e+00, 1.000000000000e+00',
       'Angle':'0.0' })
  elif Types[i]=='Z':
    MasterNode0=ET.SubElement(yarn,'MasterNode', {'index':str(ind),
       'Position': str(x0)+', '+str(y0)+', '+str(slices[i]*ImgRes),
       'Tangent':'0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00',
       'Up':'1.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00',
       'Angle':'0.0' })

  for child in yarn:
    if child.tag=='YarnSection':
       YarnSection=child
       break

PrintXML(Root,Path,OutName)
#Add nodes at the end and check direction for up vector




#Does file exist? 
#Yes:
#Does yarn exist?
#Yes:
#
#Is there 
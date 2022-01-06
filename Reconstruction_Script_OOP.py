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

#This script loads polygon xy data from ImageJ, stores and sorts them using auxiliary classes, populates TexGen classes
#and saves texgen model.

#InsertSection function will automatically initialise a new yarn for the given index and assign
#the initial section as well as position the yarn in the binary tree

#Yarns are partitioned in section sequences which determine the number of master nodes

#When inserting additional sections an interal binary tree of sections is populated
#If a different polygon is inserted for the same slice number the polygon is updated and a
#message is displayed 

def NodeDistance(N0,N1): 
  dx=N0.Position[0]-N1.Position[0]
  dy=N0.Position[1]-N1.Position[1]
  dz=N0.Position[2]-N1.Position[2]
  return m.sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2))
def ClosestIndex(Arr,TargetAngle,Centroid):
  minang=180
  minind=0
  i=0
  for point in Arr:
    v=point-Centroid
    vlen=np.linalg.norm(v)
    vu=v/vlen
    th=m.atan2(vu[0],vu[1])
    if abs(TargetAngle-th)<minang:
      minang=abs(TargetAngle-th)
      minind=i
    i+=1
  print(float(minind)/float(len(Arr)))  
  return BringToTop(Arr,minind) 

def MatchArrays(Arr1,Arr2):
  tol=1e-9
  diff_0=1.0e8
  if len(Arr1)>len(Arr2):
     diff_Arr=Arr1[:len(Arr2)]-Arr2     
  else:  
     diff_Arr=Arr1-Arr2[:len(Arr1)]       

  dist_Arr=np.sqrt(np.sum(diff_Arr**2,axis=1))
  diff_1=np.sum(dist_Arr)
  stop=0
  while tol<diff_0-diff_1 and stop<100:
    diff_0=diff_1
    Arr2=np.insert(Arr2,0,Arr2[-1],0)
    Arr2=np.delete(Arr2,-1,0)
    if len(Arr1)>len(Arr2):
       diff_Arr=Arr1[:len(Arr2)]-Arr2     
    else:  
       diff_Arr=Arr1-Arr2[:len(Arr1)]    
    dist_Arr=np.sqrt(np.sum(diff_Arr**2,axis=1))
    diff_1=np.sum(dist_Arr)   
    stop+=1
  return Arr2  
def RotateArray(Arr,Portion): # If Portion = 1  : Full rotation
  for i in range(int(len(Arr)*Portion)):
      Arr=np.insert(Arr,0,Arr[-1],0)
      Arr=np.delete(Arr,-1,0) 
  return Arr 

def BringToTop(Arr,Index): # If Portion = 1  : Full rotation
  for i in range(Index+1):
      Arr=np.insert(Arr,-1,Arr[0],0)
      Arr=np.delete(Arr,0,0) 
  return Arr 
def Centroid(Polygon):
   cx=np.sum(Polygon[:,0])/len(Polygon[:,0])
   cy=np.sum(Polygon[:,1])/len(Polygon[:,1])    
   return np.array([cx,cy])

class MasterNode:

  def __init__(self,Index,Polygon,Position,Angle,Tangent,Up):
    self.Index=Index
    self.Polygon=Polygon
    self.Position=Position
    self.Angle=Angle
    self.Tangent=Tangent
    self.Up=Up
    self.right=None
    self.left=None

  def Insert(self,Node):

    if isinstance(self.Index, int):
       
        if self.Index>Node.Index: 

           if self.left==None:
              self.left=Node
              #print str(self.Index)+' assign left node '+str(Node.Index)
           else:
              self.left.Insert(Node)

        elif self.Index<Node.Index:

           if self.right==None:
              self.right=Node
              #print str(self.Index)+' assign right node '+str(Node.Index)
           else:
              self.right.Insert(Node)

        elif self.Index==Node.Index:
           self.Position=Node.Position
           self.Angle=Node.Angle
           self.Tangent=Node.Tangent
           self.Up=Node.Up
           #print str(self.Index)+' replace node '+str(Node.Index)
    else:

        self=Node

  def GetList(self,Empty=[]):

     if self.left:
        Empty=self.left.GetList(Empty)
     
     Empty.append(self)  
     
     if self.right:
        Empty=self.right.GetList(Empty)    
     
     return Empty

  def Find(self,ind):

    if ind<self.Index:
      if self.left is None:
          return str(ind)+' Not Found'
      return self.left.Find(ind)
    elif ind>self.Index:
       if self.right is None:
          return str(ind)+' Not Found'
       return self.right.Find(ind)
    elif ind==self.Index:
       return self
  
  def CountNodes(self,Count=0):
     if self.left:
       Count=self.left.CountNodes(Count)
     Count+=1  
     if self.right:
       Count=self.right.CountNodes(Count)
     return  Count

  def PrintTree(self):
    if self.left:
      self.left.PrintTree()
    print(self.Index)
    if self.right:
      self.right.PrintTree()

if __name__=='__main__':
  #Get in the appropriate VF folder 
  cwd=cwd+'\\VF55'
  os.chdir(cwd)
  
  #Get data for window size and resolution:
  
  WinSize=np.genfromtxt('window_size.txt')
  file=open('pixel_size.txt','r')
  ImgRes=file.read() #mm
  file.close()
  ImgRes=float(ImgRes)
  #Define domain
  DS=WinSize*np.array([1.0,1.0,1.0])*ImgRes
  P0=np.array([0.0,0.0,0.0])
  CP2=XYZ(DS[0],DS[1],DS[2])
  CP1=XYZ(P0[0],P0[1],P0[2])
  CDomain=CDomainPlanes(CP1,CP2) # TexGen domain class
  #Polygon Data folder:
  DatFold='\\Data7'
  os.chdir(cwd+DatFold)
  # Store Files:
  FileList=[(f.replace('.dat','')).split('_') for f in listdir(cwd+DatFold) ] # list of info from names
  FileNames=[f for f in listdir(cwd+DatFold) ] # full names
  # Initialise auxiliary class: 
  MyYarns=Yarns(0)

  for i,file in enumerate(FileList):
    # File name structure : Direction_YarnIndex_Index_Slice ex. : X_2_1_230
    Direction=file[0]
    YarnIndex=int(file[1])
    Index=int(file[2])
    Slice=int(file[3])*ImgRes
    Sign=None
    Polygon=np.genfromtxt(cwd+DatFold+'\\'+FileNames[i])*ImgRes
    #Polygon=RotateArray(Polygon,0.5)
    if Direction in ['Z','z']:

       try: 
          Sign=int(file[4])
          Slice=DS[2]-Slice #Reversed Z axis
       except IndexError:
          #Slice=DS[2]-Slice
          pass
       Polygon=Polygon*np.array([1.0,-1.0])+np.array([0.0, DS[1]])
       C=Centroid(Polygon)
       MasterPos=np.array([C[0],C[1],Slice])

    elif Direction in ['X','x']:
       Polygon=Polygon*np.array([-1.0,-1.0])+np.array([DS[1],DS[2]])
       C=Centroid(Polygon)
       MasterPos=np.array([Slice,C[0],C[1]])

    elif Direction in ['Y','y']:
       pass         
    #Populate trees   
    MySection=Section(Index,Slice,Polygon,Direction,Sign)
    MyYarns.InsertSection(YarnIndex,MySection)

  # Add master nodes to join sections and compute final global positions  
  #MyYarns.PrintYarnTree()
  
  MyYarns.AddMasterNodes()
  
  #######
  # TexGen classes initialisation: 
  Textile=CTextile()
  Interpolation=CInterpolationBezier(False, False, False)
  
  #Traverse auxiliary yarn tree and extract yarns
  MyYarnDict=MyYarns.ExtractYarns()
  
  # Iterate yarn class to extact data and populate corresponding TexGen classes
  for y in MyYarnDict:
    MyYarn=MyYarnDict[y]
    MyNodes=MyYarn.Nodes
    NodeList=MyNodes.GetList([])
##### Just added - Add extra nodes between partition links 
    num=len(NodeList)
    if num>2:
      for i in range(num//2-1):
         N0=NodeList[3*i+1]
         N1=NodeList[3*i+2]
         dv=N1.Position-N0.Position
         NmidPos=N0.Position+dv*0.5
         tan_i=dv/np.linalg.norm(dv)
         up_i=np.cross(tan_i,np.array([0.0,1.0,0.0]))
         Nmid=MasterNode(3*i+2,NmidPos,0.0,tan_i,up_i)
         NodeList.insert(3*i+2,Nmid)
############################# 

    n0y=NodeList[0].Position[1]
    NumSlices=MyYarn.CountSlices(0)
    MySections=MyYarn.Sections
    CSection=CYarnSectionInterpPosition()
    CNodeList=[CNode(XYZ(n.Position[0],n.Position[1],n.Position[2])) for n in NodeList]
    


    for sec in MySections:
      MySection=MySections[sec]
      Direction=MySection.Direction
      index=MySection.Index
      SectionsDict=MySection.TreeToDictionary({})
      Sign=MySection.Sign
      #Local coordinate system:
      # - Adjust transformations accordingly to match the global representation
      for s in SectionsDict:
        CXYVector=XYVector()
        MyPolygon=SectionsDict[s].Polygon
        N=MyNodes.Find(2*sec)

        if Direction in ['X','x']:
          MNPos=np.array([N.Position[1],N.Position[2]])
          #MyPolygon=ClosestIndex(MyPolygon,StartAngleX,MNPos)
          LocPolygon=(MyPolygon-MNPos)*np.array([-1.0,1.0])
          LocPolygon=LocPolygon[::-1] # Fixes hollow rendering (if needed)
        elif Direction in ['Y','y']:
          MNPos=np.array([N.Position[0],N.Position[2]])
          LocPolygon=(MyPolygon-MNPos)*np.array([1.0,1.0])
          LocPolygon=LocPolygon[::-1]
        elif Direction in ['Z','z']:
          MNPos=np.array([N.Position[0],N.Position[1]])
          #MyPolygon=BringToTop(MyPolygon,3)
          if Sign:
             LocPolygon=(MyPolygon-MNPos)*np.array([1.0, -1.0])
          else:
             LocPolygon=(MyPolygon-MNPos)*np.array([-1.0, -1.0])
             LocPolygon=LocPolygon[::-1]
        else:
          print 'Unrecognised direction'   
        #LocPolygon=RotateArray(LocPolygon,0.9)   
        CXYList=[XY(p[0],p[1]) for p in LocPolygon]
        for i in CXYList:
          CXYVector.push_back(i)
        #d=SectionsDict[s].d
        CSection.AddSection(s,CSectionPolygon(CXYVector))
    CYarn0=CYarn()
    CYarn0.AssignSection(CSection)
    i=0
    for n in CNodeList:
       CYarn0.AddNode(n)
       n0=CYarn0.GetNode(i)
       Up=NodeList[i].Up
       Tangent=NodeList[i].Tangent
       CUp=XYZ(Up[0],Up[1],Up[2])        
       CTangent=XYZ(Tangent[0],Tangent[1],Tangent[2])
       try:
         n0.SetTangent(CTangent)
         n0.SetUp(CUp)
       except AttributeError:
         pass
       i+=1
    CYarn0.AssignInterpolation(Interpolation)
    if y in [33,34,35,36,37,38,39,40]:
      CYarn0.SetResolution(int(NumSlices*0.5),60)
    else:
      CYarn0.SetResolution(int(NumSlices*0.8),100) 
    Textile.AddYarn(CYarn0)
  #Save tg3 file  

  Textile.AssignDomain(CDomain)
  AddTextile('Rec4',Textile)
  SaveToXML(cwd+'\\Reconstruction5.tg3',"",OUTPUT_STANDARD)

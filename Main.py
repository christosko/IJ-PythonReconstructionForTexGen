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

class Yarns: # This class, after initialised, takes sections using InsertSection and either creates a new yarn in the tree or updates an existing one

   def __init__(self,Index):
      self.Index=Index
      self.Sections={} #Each yarn has a section list assign to it. 
      self.Nodes={} # Empty node list
      #self.Type=Type # To determine tangent direction and master node coordinates
      self.Length=0
      self.left=None
      self.right=None


   def InsertSection(self,Index,Section): #Recursive method which adds data to corresponding yarn and section 

      if isinstance(self.Index, int):
         
         if Index<self.Index:

           if self.left is None:
              self.left=Yarns(Index)              
              self.left.Sections[Section.Index]=Section
              print 'Left insert:'+str(Index)
           else:
              self.left.InsertSection(Index,Section)

         elif Index>self.Index:

           if self.right is None:
              self.right=Yarns(Index)
              self.right.Sections[Section.Index]=Section
              print 'Right insert:'+str(Index)
           else:
              self.right.InsertSection(Index,Section)

         elif Index==self.Index:

           try:
             print 'New Insert'
             self.Sections[Section.Index].Insert(Section)             
           except KeyError:
             self.Sections[Section.Index]=Section  
                   
           for S in self.Sections:

              if self.Sections[S].CountNodes(0)>1:
                  Min=self.Sections[S].FindMin()
                  Max=self.Sections[S].FindMax()
                  #print Min,Max             
                  self.Sections[S].UpdatePositions(Min.Slice,Max.Slice)

      elif self.Index==None:
         self.Index=Index
         try:
           self.Sections[Section.Index].Insert(Section)
         except KeyError:
           self.Sections[Section.Index]=Section

   def AddMasterNodes(self): # Called after data insertion is finalised

      if self.left:
        self.left.AddMasterNodes()

      for S in range(len(self.Sections)):
         MaxS=self.Sections[S].FindMax()
         MinS=self.Sections[S].FindMin()
         m0x=np.sum(MinS.Polygon[:,0])/len(MinS.Polygon[:,0])
         m0y=np.sum(MinS.Polygon[:,1])/len(MinS.Polygon[:,1])
         m1x=np.sum(MaxS.Polygon[:,0])/len(MaxS.Polygon[:,0])
         m1y=np.sum(MaxS.Polygon[:,1])/len(MaxS.Polygon[:,1])
     
         if self.Sections[S].Direction in ['X','x']:
            m0=np.array([MinS.Slice,m0x,m0y])
            m1=np.array([MaxS.Slice,m0x,m0y])
            t=np.array([1.0,0.0,0.0])
            up=np.array([0.0,0.0,1.0])
            M0=MasterNode(0,m0,0.0,t,up)
            M1=MasterNode(1,m1,0.0,t,up)    
         elif self.Sections[S].Direction in ['Y','y']:
            m0=np.array([m0x,MinS.Slice,m0y])
            m1=np.array([m0x,MaxS.Slice,m0y])
            t=np.array([0.0,1.0,0.0])
            up=np.array([0.0,0.0,1.0])
            M0=MasterNode(0,m0,0.0,t,up)
            M1=MasterNode(1,m1,0.0,t,up)
         elif self.Sections[S].Direction in ['Z','z']:
            m1=np.array([m1x,m1y,MinS.Slice])
            m0=np.array([m1x,m1y,MaxS.Slice])
            t=np.array([0.0,0.0,1.0])
            up=np.array([1.0,0.0,0.0])
            M0=MasterNode(0,m0,0.0,t,up)
            M1=MasterNode(1,m1,0.0,t,up)     
         else :
            print 'No direction specified for Section:'+str(self.Sections[S].Index)     

         self.Nodes[2*S]=M0
         self.Nodes[2*S+1]=M1

      for i in range(len(self.Nodes)-1):
         self.Length+=NodeDistance(self.Nodes[i],self.Nodes[i+1])

      print self.Length
      if len(self.Sections)>1:
         PreLength=0 
         for S in range(len(self.Sections)):
           SecLength=NodeDistance(self.Nodes[2*S],self.Nodes[2*S+1])
           self.Sections[S].UpdateGlobalPositions(self.Length,PreLength,SecLength)
           try:
             PreLength+=SecLength+NodeDistance(self.Nodes[2*S+1],self.Nodes[2*S+2])
           except KeyError:
             pass  

      if self.right:
        self.right.AddMasterNodes()

   def ExtractSections(self,Index): 

     if isinstance(self.Index, int):
       if Index<self.Index:
          if self.left is None:
             return 'Yarn'+str(Index)+'Not Found'
          return self.left.ExtractSection(Index)
       elif Index>self.Index:
          if self.right is None:
             return 'Yarn'+str(Index)+'Not Found'
          return self.right.ExtractSection(Index)
       else:
          return self.Sections

   def PrintYarnTree(self):

      if self.left:
        self.left.PrintYarnTree()
      print(self.Index)
      if self.right:
        self.right.PrintYarnTree()

   def ExtractYarns(self,EmptyDict):
     if self.left:
       self.left.ExtractYarns(EmptyDict)
     EmptyDict[self.Index]=self  
     if self.right:
       self.right.ExtractYarns(EmptyDict)    
     return EmptyDict

class Section:

  def __init__(self,Index,Slice,Polygon,Direction):

    self.Index=Index
    self.Slice=Slice
    self.Direction=Direction
    self.d=0.0
    self.Polygon=Polygon
    if self.Slice:
       CentrePos=np.sum(Polygon,0)/len(Polygon[:,0])
       if Direction in ['X','x']:
          Position=np.array([Slice,CentrePos[0],CentrePos[1]])
          self.SlaveNode=SlaveNode(Position)
       elif Direction in ['Y','y']:
          Position=np.array([CentrePos[0],Slice,CentrePos[1]])
          self.SlaveNode=SlaveNode(Position)        
       elif Direction in ['Z','z']:   
          Position=np.array([CentrePos[0],CentrePos[1],Slice])
          self.SlaveNode=SlaveNode(Position)   
    else:
       self.SlaveNode=None     
       
    self.left=None
    self.right=None
    
  def Insert(self,Section): #New node in section tree

    if isinstance(self.Slice,float or int):
         if self.Slice>Section.Slice: 

           if self.left==None:
              print 'Left section insert:'+str(Section.Slice)
              self.left=Section
           else:
              self.left.Insert(Section)

         elif self.Slice<Section.Slice:

           if self.right==None:
              print 'Right section insert:'+str(Section.Slice)
              self.right=Section
           else:
              self.right.Insert(Section)

         elif self.Slice==Section.Slice:
           print('Polygon updated for slice:'+str(Section.Slice)+' in Yarn:'+str(Section.Index))
           self.Polygon=Section.Polygon

    elif self.Slice==None:  
      self=Section


  def UpdatePositions(self,Min,Max):

    if isinstance(self.Slice,float or int):
      if self.left:
        self.left.UpdatePositions(Min,Max) 
      self.d=float((self.Slice-Min))/float((Max-Min))
      if self.right:
        self.right.UpdatePositions(Min,Max)

  def UpdateGlobalPositions(self,YarnLength,PreLength,SecLength): # For entire path
   
   if self.left:
     self.left.UpdateGlobalPositions(YarnLength,PreLength,SecLength) 
   self.d=(PreLength+SecLength*self.d)/YarnLength
   if self.right:
     self.right.UpdateGlobalPositions(YarnLength,PreLength,SecLength)    
 
  def FindMax(self):
   if isinstance(self.Slice, float or int):
       if self.right: 
          return self.right.FindMax()
       return self

  def FindMin(self):
   if isinstance(self.Slice,float or int):
       if self.left: 
          return self.left.FindMin()
       return self

  def CountNodes(self,Count):
     if self.left:
       Count=self.left.CountNodes(Count)
     Count+=1  
     if self.right:
       Count=self.right.CountNodes(Count)
     return  Count
    
  def TreeToDictionary(self,EmptyDict):
     if self.left:
       EmptyDict=self.left.TreeToDictionary(EmptyDict)

     EmptyDict[self.d]=self  

     if self.right:
       EmptyDict=self.right.TreeToDictionary(EmptyDict)
     return  EmptyDict

  def PrintTree(self):
     if self.left:
       self.left.PrintTree()
     print(self.Slice)
     if self.right:
       self.right.PrintTree()  

class MasterNode:
  def __init__(self,Index,Position,Angle,Tangent,Up):
    self.Index=Index
    self.Position=Position
    self.Angle=Angle
    self.Tangent=Tangent
    self.Up=Up
 
class SlaveNode:
  def __init__(self,Position):
    self.Position=Position


if __name__=='__main__':
  #####
  
  #Get data for window size and resolution:
  WinSize=np.genfromtxt('window_size.txt')
  file=open('pixel_size.txt','r')
  ImgRes=file.read() #mm
  file.close()
  ImgRes=float(ImgRes)
  #Define domain
  DS=WinSize*np.array([1.0,1.0,-1.0])*ImgRes
  print DS
  P0=np.array([0.0,0.0,0.0])
  CP2=XYZ(DS[0],DS[1],DS[2])
  CP1=XYZ(P0[0],P0[1],P0[2])
  CDomain=CDomainPlanes(CP1,CP2) # TexGen domain class
  #Polygon Data folder:
  DatFold='\\Data'
  os.chdir(cwd+DatFold)
  # Store Files:
  FileList=[(f.replace('.dat','')).split('_') for f in listdir(cwd+DatFold)] # list of info from names
  FileNames=[f.replace('.dat','.dat') for f in listdir(cwd+DatFold)] # full names
  # Initialise auxiliary class: 
  MyYarns=Yarns(0)
  i=0
  for file in FileList:
    # File name structure : Direction_YarnIndex_Index_Slice ex. : X_2_1_230
    Direction=file[0]
    YarnIndex=int(file[1])
    Index=int(file[2])
    if Direction in ['Z','z']:
       Slice=-int(file[3])*ImgRes
       Polygon=np.genfromtxt(cwd+DatFold+'\\'+FileNames[i])*ImgRes
    else:
       Slice=int(file[3])*ImgRes
       Polygon=np.genfromtxt(cwd+DatFold+'\\'+FileNames[i])*np.array([1.0,-1.0])*ImgRes   
    #Populate trees   
    MySection=Section(Index,Slice,Polygon,Direction)
    MyYarns.InsertSection(YarnIndex,MySection)
    i+=1
  # Add master nodes to join sections and compute final global positions  
  #MyYarns.PrintYarnTree()
  
  MyYarns.AddMasterNodes()

  # TexGen classes initialisation: 
  Textile=CTextile()
  Interpolation=CInterpolationLinear(False, False, False)
  
  #Traverse auxiliary yarn tree and extract yarns
  MyYarnDict=MyYarns.ExtractYarns({})
  # Iterate yarn class to extact data and populate corresponding TexGen classes
  for y in MyYarnDict:
    MyYarn=MyYarnDict[y]
    Nodes=MyYarn.Nodes
    print Nodes
    MySections=MyYarn.Sections
    CSection=CYarnSectionInterpPosition()
    CNodeList=[CNode(XYZ(Nodes[n].Position[0],Nodes[n].Position[1],Nodes[n].Position[2])) for n in Nodes]
    
    for sec in MySections:
      MySection=MySections[sec]
      Direction=MySection.Direction
      index=MySection.Index
      SectionsDict=MySection.TreeToDictionary({})
      for s in SectionsDict:
        CXYVector=XYVector()
        MyPolygon=SectionsDict[s].Polygon
        print MyPolygon
        if Direction in ['X','x']:
          MNPos=np.array([Nodes[2*sec].Position[1],Nodes[2*sec].Position[2]])
          print MNPos
          LocPolygon=MyPolygon-MNPos
        elif Direction in ['Y','y']:
          MNPos=np.array([Nodes[2*sec].Position[0],Nodes[2*sec].Position[2]])
          print MNPos
          LocPolygon=MyPolygon-MNPos
        elif Direction in ['Z','z']:
          MNPos=np.array([Nodes[2*sec].Position[0],Nodes[2*sec].Position[1]])
          print MNPos
          LocPolygon=MyPolygon-MNPos
        else:
          print 'Unrecognised direction'    
        CXYList=[XY(p[0],p[1]) for p in LocPolygon]
        for i in CXYList:
          CXYVector.push_back(i)
        d=SectionsDict[s].d
        CSection.AddSection(d,CSectionPolygon(CXYVector))
    CYarn0=CYarn()
    CYarn0.AssignSection(CSection)
    i=0
    for N in CNodeList:
       CYarn0.AddNode(N)
       n0=CYarn0.GetNode(i)
       Up=Nodes[i].Up
       CUp=XYZ(Up[0],Up[1],Up[2])
       n0.SetUp(CUp)
       i+=1
    CYarn0.AssignInterpolation(Interpolation)
    CYarn0.SetResolution(50,100)
    Textile.AddYarn(CYarn0)
  #Save tg3 file  
  Textile.AssignDomain(CDomain)
  AddTextile('Rec',Textile)
  SaveToXML(cwd+'\\Reconstruction.tg3',"",OUTPUT_STANDARD)

    

    




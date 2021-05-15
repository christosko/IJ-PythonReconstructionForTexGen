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

class Yarns: # This class, after initialised, takes sections using InsertSection and either creates a new yarn in the tree or updates an existing one

   def __init__(self,Index):
      self.Index=Index
      self.Sections={} #Each yarn has a section list assigned to it. 
      self.Nodes=MasterNode(0,None,None,None,None) # Empty node class
      self.Length=0
      self.left=None
      self.right=None


   def InsertSection(self,Index,Section): #Recursive method which adds data to corresponding yarn and section 

      if isinstance(self.Index, int):
         
         if Index<self.Index:

           if self.left is None:
              self.left=Yarns(Index)              
              self.left.Sections[Section.Index]=Section
              #print 'Left insert:'+str(Index)
           else:
              self.left.InsertSection(Index,Section)

         elif Index>self.Index:

           if self.right is None:
              self.right=Yarns(Index)
              self.right.Sections[Section.Index]=Section
              #print 'Right insert:'+str(Index)
           else:
              self.right.InsertSection(Index,Section)

         elif Index==self.Index:

           try:
             #print 'New Insert'
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

      for S in self.Sections:
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
            m0=np.array([m0x,m0y,MinS.Slice])
            m1=np.array([m0x,m0y,MaxS.Slice])
            if self.Sections[S].Sign:
               t=np.array([0.0,0.0,1.0])              
               up=np.array([0.0,-1.0,0.0])                 
               M0=MasterNode(0,m0,0.0,t,up)
               M1=MasterNode(1,m1,0.0,t,up) 
            else:
               t=np.array([0.0,0.0,-1.0])              
               up=np.array([0.0,-1.0,0.0])              
               M0=MasterNode(0,m1,0.0,t,up)
               M1=MasterNode(1,m0,0.0,t,up)                                 
         else :
            print 'No direction specified for Section:'+str(self.Sections[S].Index)  
         # Rotates polygon point order to match the previous polygon -  not great
         #if S==0:
         #   self.Sections[S].MatchPolygons(MinS.Polygon)
         #else:
         #   self.Sections[S].MatchPolygons(self.Sections[S-1].FindMax().Polygon)
         M0.Index=2*S
         M1.Index=2*S+1 
        #print(M0.Position, M1.Position)
         self.Nodes.Insert(M0)
         self.Nodes.Insert(M1)
         #print 'In yarn: '+str(self.Index)+' and section: '+str(S)+' nodes: '+str(2*S)+','+str(2*S+1)
      # Fix duplicate node issue
      Nodes=self.Nodes
      NodeList=Nodes.GetList([])

      for i in range(Nodes.CountNodes()-1):
         self.Length+=NodeDistance(NodeList[i],NodeList[i+1])
      SecInd=sorted([i for i in self.Sections])
      #print self.Length
      if len(self.Sections)>1:
         PreLength=0 
         for i in range(len(self.Sections)):
           print(NodeList[2*i].Index,NodeList[2*i+1].Index)
           SecLength=NodeDistance(NodeList[2*i],NodeList[2*i+1])
           self.Sections[SecInd[i]].UpdateGlobalPositions(self.Length,PreLength,SecLength)
           try:
             PreLength+=SecLength+NodeDistance(NodeList[2*i+1],NodeList[2*i+2])
           except IndexError:
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

   def ExtractYarns(self,EmptyDict={}):
     if self.left:
       self.left.ExtractYarns(EmptyDict)
     EmptyDict[self.Index]=self  
     if self.right:
       self.right.ExtractYarns(EmptyDict)    
     return EmptyDict

   def CountSlices(self,Count=0):

     if self.left:
        self.left.CountSlices(Count)
 
     for sec in self.Sections:
        Count=self.Sections[sec].CountNodes(Count)
 
     if self.right:
        self.right.CountSlices(Count)
     
     return Count
 
       

class Section: #Data binary tree for signle direction section sequence - recursive method

  def __init__(self,Index,Slice,Polygon,Direction,Sign):

    self.Index=Index
    self.Slice=Slice
    self.Direction=Direction
    self.d=0.0
    self.Polygon=Polygon
    self.Sign=Sign
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
              self.left=Section
           else:
              self.left.Insert(Section)

         elif self.Slice<Section.Slice:

           if self.right==None:
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
  
  def MatchPolygons(self,PreviousPolygon): 
    if self.left:
       self.left.MatchPolygons(self.Polygon)
    self.Polygon=MatchArrays(PreviousPolygon,self.Polygon) 
    if self.right:
       self.right.MatchPolygons(self.Polygon)    
 
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

  def CountNodes(self,Count=0):
     if self.left:
       Count=self.left.CountNodes(Count)
     Count+=1  
     if self.right:
       Count=self.right.CountNodes(Count)
     return  Count
    
  def TreeToDictionary(self,EmptyDict={}):
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

class SlaveNode:
  def __init__(self,Position):
    self.Position=Position


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
  DatFold='\\Data3'
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
    Slice=int(file[3])*ImgRes
    Sign=None
    Polygon=np.genfromtxt(cwd+DatFold+'\\'+FileNames[i])*ImgRes
    #Polygon=RotateArray(Polygon,0.5)
    if Direction in ['Z','z']:
       try: 
          Sign=int(file[4])
          Slice=DS[2]-Slice
       except IndexError:
          #Slice=DS[2]-Slice
          pass
       Polygon=Polygon*np.array([1.0,-1.0])+np.array([0.0,DS[1]])#([DS[0],DS[1]])
    elif Direction in ['X','x']:
       #Slice=DS[0]-Slice
       Polygon=Polygon*np.array([-1.0,-1.0])+np.array([DS[1],DS[2]])
    elif Direction in ['Y','y']:
       Polygon=Polygon*np.array([1.0,-1.0])+np.array([0.0,DS[2]])             
    #Populate trees   
    MySection=Section(Index,Slice,Polygon,Direction,Sign)
    MyYarns.InsertSection(YarnIndex,MySection)
    i+=1
  # Add master nodes to join sections and compute final global positions  
  #MyYarns.PrintYarnTree()
  
  MyYarns.AddMasterNodes()

  # TexGen classes initialisation: 
  Textile=CTextile()
  Interpolation=CInterpolationLinear(False, False, False)
  
  #Traverse auxiliary yarn tree and extract yarns
  MyYarnDict=MyYarns.ExtractYarns()
  # Iterate yarn class to extact data and populate corresponding TexGen classes
  for y in MyYarnDict:
    MyYarn=MyYarnDict[y]
    MyNodes=MyYarn.Nodes
    NodeList=MyNodes.GetList([])
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
          LocPolygon=(MyPolygon-MNPos)*np.array([-1.0,1.0])
          LocPolygon=LocPolygon[::-1] # Fixes hollow rendering (if needed)
        elif Direction in ['Y','y']:
          MNPos=np.array([N.Position[0],N.Position[2]])
          LocPolygon=(MyPolygon-MNPos)*np.array([1.0,1.0])
          LocPolygon=LocPolygon[::-1]
        elif Direction in ['Z','z']:
          MNPos=np.array([N.Position[0],N.Position[1]])
          if Sign:
             LocPolygon=(MyPolygon-MNPos)*np.array([1.0,-1.0])
          else:
             LocPolygon=(MyPolygon-MNPos)*np.array([-1.0,-1.0])
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
    CYarn0.SetResolution(int(NumSlices*0.6),100)
    Textile.AddYarn(CYarn0)
  #Save tg3 file  

  Textile.AssignDomain(CDomain)
  AddTextile('Rec',Textile)
  SaveToXML(cwd+'\\Reconstruction2.tg3',"",OUTPUT_STANDARD)


    





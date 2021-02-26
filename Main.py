import numpy as np

import math as m
import os

cwd=os.getcwd()
sys.path.append(cwd)
#Requires a minimum build of texgen and a path pointing to Python libraris
TexGenPath='C:\\Python27\\Lib\\site-packages\\TexGen'
sys.path.append(TexGenPath)

from TexGen.Core import*

#Auxiliary class is used to sort yarns in binary tree, assign attributes and sections
#InsertSection function will automatically initialise a new yarn for the given index and assign
#the initial section as well as position the yarn in the binary tree
#Section class is initialised with a slice number and polygon array which are obtained from each corresponding 
#dat. file 
#When inserting additional sections an interal binary tree of sections is populated
#If a different polygon is inserted for the same slice number the polygon is updated and a
#message is displayed 

def NodeDistance(N0,N1): 
  dx=N0.Position[0]-N1.Position[0]
  dy=N0.Position[1]-N1.Position[1]
  dz=N0.Position[2]-N1.Position[2]
  return sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2))
class Yarns: # This class, after initialised, takes sections using InsertSection and either creates a new yarn in the tree or updates an existing one

   def __init__(self,Index):
      self.Index=Index
      self.Sections={0:Section(0,None,None,None)} #Each yarn has a section list assign to it. 
      self.Nodes={} # Empty node list
      #self.Type=Type # To determine tangent direction and master node coordinates
      self.Length=0
      self.left=None
      self.right=None


   def InsertSection(self,Index,Section): #Recursive method which adds (Slice,Polygon) data to corresponding yarn and section 

      if isinstance(self.Index,int):
         
         if Index<self.Index:

           if self.left is None:
              self.left=Yarns(Index)              
              try:
                self.left.Sections[Section.Index].Insert(Section)     
              except KeyError:
                self.left.Sections[Section.Index]=Section
           else:
              self.left.InsertSection(Index,Section)

         elif Index>self.Index:

           if self.right is None:
              self.right=Yarns(Index)
              try:
                self.right.Sections[Section.Index].Insert(Section)
              except KeyError:
                self.Sections[Section.Index]=Section
           else:
              self.right.InsertSection(Index,Section)

         elif Index==self.Index:
           if self.Section.Slice:
              try:
                self.Sections[Section.Index].Insert(Section)             
              except KeyError:
                self.Sections[Section.Index]=Section  
              USections=self.ExtractSections(Index)
              for S in USections:
                Min=USections[S].FindMin()
                Max=USections[S].FindMax()
                #print Min,Max
                self.Sections[S].UpdatePositions(Min,Max)
              self.AddNodes()
           else:
              self.Sections[Section.Index]=Section           
      elif self.Index==None:
         self.Index=Index
         try:
           self.Sections[Section.Index].Insert(Section)
         except KeyError:
           self.Sections[Section.Index]=Section

   def AddMasterNodes(self):

      if self.left:
        self.left.AddMasterNodes()
      for S in range(len(self.Sections)):
        MaxS=Sections[S].FindMax()
        MinS=Sections[S].FindMin()
        m0x=np.sum(MinS.Polygon[:,0])/len(MinS.Polygon[:,0])
        m0y=np.sum(MinS.Polygon[:,1])/len(MinS.Polygon[:,1])
        m1x=np.sum(MaxS.Polygon[:,0])/len(MaxS.Polygon[:,0])
        m1y=np.sum(MaxS.Polygon[:,1])/len(MaxS.Polygon[:,1])
    
        if Sections[S].Direction in ['X','x']:
           m0=np.array([MinS.Slice,m0x,m0y])
           m1=np.array([MaxS.Slice,m0x,m0y])
           t=np.array([1.0,0.0,0.0])
           up=np.array([0.0,0.0,1.0])
           M0=MasterNode(0,m0,0.0,t,up)
           M1=MasterNode(1,m1,0.0,t,up)    
        elif Sections[S].Direction in ['Y','y']:
           m0=np.array([m0x,MinS.Slice,m0y])
           m1=np.array([m0x,MaxS.Slice,m0y])
           t=np.array([0.0,1.0,0.0])
           up=np.array([0.0,0.0,1.0])
           M0=MasterNode(0,m0,0.0,t,up)
           M1=MasterNode(1,m1,0.0,t,up)
        elif Sections[S].Direction in ['Z','z']:
           m0=np.array([m0x,m0y,MinS.Slice])
           m1=np.array([m0x,m0y,MaxS.Slice])
           t=np.array([0.0,0.0,1.0])
           up=np.array([1.0,0.0,0.0])
           M0=MasterNode(0,m0,0.0,t,up)
           M1=MasterNode(1,m1,0.0,t,up)            
        self.Nodes[2*S]=M0
        self.Nodes[2*S+1]=M1
           
      FullLength=0
      for i in range(len(self.Nodes)-1):
         FullLength+=NodeDistance(self.Nodes[i],self.Nodes[i+1])
      self.Length=FullLength  
      PreLength=0 
      for S in range(len(self.Sections)):
        SecLength=NodeDistance(self.Nodes[2*S],self.Nodes[2*S+1])
        self.Sections[S].UpdateGlobalPositions(FullLength,PreLength,SecLength)
        PreLength+=SecLength+NodeDistance(self.Nodes[2*S+1],self.Nodes[2*S+2])

      if self.right:
        self.right.AddMasterNodes()

   def ExtractSections(self,Index):

     if isinstance(self.Index,int):
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
      print(self.Section)
      if self.right:
        self.right.PrintYarnTree()


   def ExtractYarns(self,EmptyDict):
     if self.left:
       self.left.ExtractYarns()
     EmptyDict[self.Index]=self
     if self.right:
       self.right.ExtractYarns()    


class Section:

  def __init__(self,Index,Slice,Polygon,Direction):
    self.Index=Index
    self.Slice=Slice
    self.Direction=Direction
    self.d=0.0
    self.Polygon=Polygon
    if self.Polygon:
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
    
  def Insert(self,Section):
    if self.Index==Section.Index:
      if isinstance(self.Slice,int):
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

    if isinstance(self.Slice,int):
      if self.left:
        self.left.UpdatePositions(Min,Max) 
      self.d=float((self.Slice-Min))/float((Max-Min))
      if self.right:
        self.right.UpdatePositions(Min,Max)

  def UpdateGlobalPositions(self,YarnLength,PreLength,SecLength):
   
   if self.left:
     self.left.UpdatePositions(Min,Max) 
   self.d=(PreLength+SecLength*self.t)/YarnLength
   if self.right:
     self.right.UpdatePositions(Min,Max)    
 
  def FindMax(self):
    if isinstance(self.Slice,int):
       if self.right: 
          return self.right.FindMax()
       return self.Slice

  def FindMin(self):
    if isinstance(self.Slice,int):
       if self.left: 
          return self.left.FindMin()
       return self.Slice  

  def ExtractData(self,DataDict):
     if isinstance(self.Slice,int):
        if self.left:
           DataDict=self.left.ExtractData(DataDict)
        DataDict[self.t]=self.Polygon
        if self.right:
           DataDict=self.right.ExtractData(DataDict)
     return DataDict   

  def TreeToDictionary(self,EmptyDict):
     if self.left:
       self.left.TreeToDictionary()
     EmptyDict[self.Index]=self
     if self.right:
       self.right.TreeToDictionary()

class MasterNode:
  def __init__(self,Index,Position,Angle,Tangent,Up):
    self.Index=Index
    self.Position=Position
    self.Angle=Angle
    self.Tangent=Tangent
    self.Up=Up
 
class SlaveNode:
  def __init__(self,Position)
    self.Position=Position


if __name__=='__main__':
  #####
 # TPolygon=np.array([[3.4,5],[5,6],[4.5,3.2]])
 # S0=234
 # S1=345
 # S2=400
 # S3=534
 # Y=Yarns(0)
 # Sec1=Section(S0,TPolygon)
 # Sec2=Section(S1,TPolygon)
 # Sec3=Section(S2,TPolygon)
 # Sec4=Section(S3,TPolygon)
 # Y.InsertSection(0,Sec1)
 # Y.InsertSection(0,Sec2)
 # Y.InsertSection(0,Sec3)
 # Y.InsertSection(0,Sec4)
 # Y.InsertSection(1,Sec1)
 # Y.InsertSection(1,Sec2)
 # Y.InsertSection(1,Sec3)
 # Y.InsertSection(1,Sec4)
 # Y.InsertSection(2,Sec1)
 # Y.InsertSection(2,Sec2)
 # Y.InsertSection(2,Sec3) 
 # Y.InsertSection(2,Sec4)
 # Y.PrintYarnTree()
 # Sec=Y.ExtractSection(1)
 # Sec.PrintSectionTree()
 # Dat={}
 # Dat=Sec.ExtractData(Dat)
 # print Dat
  
  #Get data for window size and resolution:
  WinSize=np.genfromtxt('window_size.txt')
  file=open('pixel_size.txt','r')
  ImgRes=file.read() #mm
  file.close()
  ImgRes=float(ImgRes)
  #Define domain
  DS=WinSize*np.array([1.0,1.0,-1.0])*ImgRes
  P0=np.array([0.0,0.0,0.0])
  CDomain=CDomainPlanes(P0,DS) # TexGen domain class
  #Polygon Data folder:
  DatFold='\\Data'
  # Store Files:
  FileList=[(f.replace('.dat','')).split('_') for f in listdir(cwd+DatFold)] # list of info from names
  FileNames=[f for f in listdir(cwd+DatFold)] # full names
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
       Polygon=np.genfromtxt(FileNames[i])*ImgRes
    else:
       Slice=int(file[3])*ImgRes
       Polygon=np.genfromtxt(FileNames[i])*np.array([1.0,-1.0])*ImgRes   
    #Populate trees   
    Section=Section(Index,Slice,Polygon,Direction)
    MyYarns.InsertSection(YarnIndex,Section)
    i+=1
  # Add master nodes to join sections and compute final global positions  
  MyYarns.AddMasterNodes()

  # TexGen classes initialisation: 
  Textile=CTextile()
  Interpolation=CInterpolationLinear(False, False, False)
  
  #Traverse auxiliary yarn tree and extract yarns
  YarnDict={}
  MyYarnDict=MyYarns.ExtractYarns(YarnDict)
  # Iterate yarn class to extact data and populate corresponding TexGen classes
  for y in MyYarnDict:
    MyYarn=MyYarnDict[y]
    Nodes=MyYarn.Nodes
    MySections=MyYarn.Sections
    CSection=CYarnSectionInterpPosition()
    CNodeList=[CNode(XYZ(Nodes[n][0],Nodes[n][1],Nodes[n][2])) for n in Nodes]
    for sec in MySections:
      MySection=MySections[sec]
      direction=MySection.Direction
      index=MySection.Index
      EmptyDict={}
      SectionsDict=MySection.TreeToDictionary(EmptyDict)
      for s in SectionsDict:
        CXYVector=XYVector()
        Polygon=SectionsDict[s].Polygon
        if Direction in ['X','x']:
          LocPolygon=Polygon-np.array([Nodes[2*sec].Position[1],Nodes[2*sec].Position[2]])
        elif Direction in ['Y','y']:
          LocPolygon=Polygon-np.array([Nodes[2*sec].Position[0],Nodes[2*sec].Position[2]])
        elif Direction in ['Z','z']:
          LocPolygon=Polygon-np.array([Nodes[2*sec].Position[0],Nodes[2*sec].Position[1]])
        else:
          print 'Unrecognised direction'  
        CXYList=[XY(p[0],p[1]) for p in LocPolygon]
        for i in CXYList:
          CXYVector.push_back(i)
        d=SectionsDict[s].d
        CSection.AddSection(d,CSectionPolygon(CXYVector))
    CYarn=CYarn()
    CYarn.AssignSection(CSection)
    for N in CNodeList:
       CYarn.AddNode(N)
    CYarn.AssignInterpolation(Interpolation)
    CYarn.SetResolution(50,100)
    Textile.AddYarn(CYarn)
  Textile.AssignDomain(CDomain)
  AddTextile('Rec',Textile)
  SaveToXML(cwd+'\\Reconstruction.tg3',"",OUTPUT_STANDARD)

    

    




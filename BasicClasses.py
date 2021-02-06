import numpy as np

import math
import os

class Yarns:

   def __init__(self,Index,Type):
      self.Index=Index
      self.Section=Section(None,None)
      self.Node=MasterNode(None,None,None,None)
      self.Type=Type # This will be used to determine master node definition method
      self.left=None
      self.right=None
      
   def InsertSection(self,Index,Section):

      if isinstance(self.Index,int):
         
         if Index<self.Index:

           if self.left is None:
             self.left=Yarns(Index)              
             self.left.Section.Insert(Index, Section.Slice, Section.Polygon)
           else:
             self.left.InsertSection(Index,Section)

         elif Index>self.Index:

           if self.right is None:
             self.right=Yarns(Index)
             self.right.Section.Insert(Index, Section.Slice, Section.Polygon)
           else:
             self.right.InsertSection(Index,Section)

         elif Index==self.Index:
           if self.Section.Slice:
              self.Section.Insert(Index, Section.Slice, Section.Polygon)             
              USec=self.ExtractSection(Index)
              Min=USec.FindMin()
              Max=USec.FindMax()
              print Min,Max
              self.Section.UpdatePositions(Min,Max)
           else:
              self.Section=Section            
      elif self.Index==None:
         self.Index=Index
         self.Section.Insert(Index, Section.Slice, Section.Polygon)     

   def AddNodes(self,PositionList):
    if PositionList:
      if self.Type=='X':
      elif self.Type=='Y':
       
         
     




   def ExtractSection(self,Index):

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
         return self.Section

   def PrintYarnTree(self):

      if self.left:
        self.left.PrintYarnTree()
      print(self.Section)
      if self.right:
        self.right.PrintYarnTree()

class Section:

  def __init__(self,Slice,Polygon):
    self.Index=None
    self.Slice=Slice
    self.t=0.0
    self.Polygon=Polygon
    self.left=None
    self.right=None
    
  def Insert(self,Index,Slice,Polygon):

    if isinstance(self.Slice,int):
      if self.Slice>Slice:
        if self.left==None:
          self.left=Section(Slice,Polygon)
          self.left.Index=Index
        else:
          self.left.Insert(Index, Slice,Polygon)

      elif self.Slice<Slice:

        if self.right==None:
          self.right=Section(Slice, Polygon)
          self.right.Index=Index
        else:
          self.right.Insert(Index, Slice, Polygon)

      elif self.Slice==Slice:
        print('Polygon updated for slice:'+str(Slice)+' in Yarn:'+str(Index))
        self.Index=Index
        self.Slice=Slice
        self.Polygon=Polygon

    elif self.Slice==None:  
      self.Index=Index
      self.Slice=Slice
      self.Polygon=Polygon

  def UpdatePositions(self,Min,Max):

    if isinstance(self.Slice,int):
      if self.left:
        self.left.UpdatePositions(Min,Max) 
      self.t=float((self.Slice-Min))/float((Max-Min))
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

  def PrintSectionTree(self):
     if self.left:
       self.left.PrintSectionTree()
     print(self.t)
     if self.right:
       self.right.PrintSectionTree()

class MasterNode:
  def __init__(self,Index,Position,Angle,Tangent,Up):
    self.Index=Index
    self.Position=Position
    self.Angle=Angle
    self.Tangent=Tangent
    self.Up=Up
    
if __name__=='__main__':
  #####
  TPolygon=np.array([[3.4,5],[5,6],[4.5,3.2]])
  S0=234
  S1=345
  S2=400
  S3=534
  Y=Yarns(0)
  Sec1=Section(S0,TPolygon)
  Sec2=Section(S1,TPolygon)
  Sec3=Section(S2,TPolygon)
  Sec4=Section(S3,TPolygon)
  Y.InsertSection(0,Sec1)
  Y.InsertSection(0,Sec2)
  Y.InsertSection(0,Sec3)
  Y.InsertSection(0,Sec4)
  Y.InsertSection(1,Sec1)
  Y.InsertSection(1,Sec2)
  Y.InsertSection(1,Sec3)
  Y.InsertSection(1,Sec4)
  Y.InsertSection(2,Sec1)
  Y.InsertSection(2,Sec2)
  Y.InsertSection(2,Sec3) 
  Y.InsertSection(2,Sec4)
  Y.PrintYarnTree()
  Sec=Y.ExtractSection(1)
  Sec.PrintSectionTree()
  Dat={}
  Dat=Sec.ExtractData(Dat)
  print Dat
# New reconstruction scheme : 
# Trace (point and click) binder yarn path in mid xz plane: 
# Information obtained: X and Z slice pairs and orientation vectors
#1) Save X,Z pairs
#2) Compute orientation vectors
#3) Save .raw files with corresponding slices for X and Z
#4) Trace 
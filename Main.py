from locale import normalize
import numpy as np
import math as m
import os
import sys
from os import listdir
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
def CheckX(V):
  Xunit=XYZ(1.0,0.0,0.0)
  Yunit=XYZ(0.0,1.0,0.0)
  Zunit=XYZ(0.0,0.0,1.0)
  if GetLength(CrossProduct(Xunit,V))<GetLength(CrossProduct(Zunit,V)) and GetLength(CrossProduct(Xunit,V))<GetLength(CrossProduct(Yunit,V)):
     return True
  else:
     return False

def CheckY(V):
  Xunit=XYZ(1.0,0.0,0.0)
  Yunit=XYZ(0.0,1.0,0.0)
  Zunit=XYZ(0.0,0.0,1.0)
  if GetLength(CrossProduct(Yunit,V))<GetLength(CrossProduct(Zunit,V)) and GetLength(CrossProduct(Yunit,V))<GetLength(CrossProduct(Xunit,V)):
     return True
  else:
     return False

def CheckZ(V):
  Xunit=XYZ(1.0,0.0,0.0)
  Yunit=XYZ(0.0,1.0,0.0)
  Zunit=XYZ(0.0,0.0,1.0)
  if GetLength(CrossProduct(Zunit,V))<GetLength(CrossProduct(Xunit,V)) and GetLength(CrossProduct(Zunit,V))<GetLength(CrossProduct(Yunit,V)):
     return True
  else:
     return False

def Check2DX(V):
  Xunit=XYZ(1.0,0.0,0.0)
  Yunit=XYZ(0.0,1.0,0.0)
  if GetLength(CrossProduct(Xunit,V))<GetLength(CrossProduct(Yunit,V)):
     return True
  else:
     return False

def Check2DY(V):
  Xunit=XYZ(1.0,0.0,0.0)
  Yunit=XYZ(0.0,1.0,0.0)
  if GetLength(CrossProduct(Yunit,V))<GetLength(CrossProduct(Xunit,V)):
     return True
  else:
     return False

def TolCheckX(V,tolerance):
   Xunit=XYZ(1.0,0.0,0.0)
   Zunit=XYZ(0.0,0.0,1.0)   
   if abs(V.z)<abs(V.x)*tolerance:
   #if GetLength(CrossProduct(Xunit,V))<GetLength(CrossProduct(Zunit,V))*tolerance:
      return True
   else:
      return False   
def TolCheckY(V,tolerance):
   Yunit=XYZ(0.0,1.0,0.0)
   Zunit=XYZ(0.0,0.0,1.0)   
   if abs(V.z)<abs(V.y)*tolerance:
   #if GetLength(CrossProduct(Yunit,V))<GetLength(CrossProduct(Zunit,V))*tolerance:
      return True
   else:
      return False  

def GetX(V):
   return V[0].x
def GetY(V):
   return V[0].y

def GetCentroidsXYZ(Yarn):
   CentVec=XYZVector()
   SlaveNodes=Yarn.GetSlaveNodes(2)
   for s in SlaveNodes:
      Orig=XYZ(0.0,0.0,0.0)
      secpts=s.GetSectionPoints()
      for p in secpts:
         Orig=Orig+p
      Orig=Orig*float(1.0/float(len(secpts)))
      CentVec.push_back(Orig)
   return CentVec   

def YarnTypeSort(Textile):
   Yarns=Textile.GetYarns()
   X=[]
   Y=[]
   Z=[]
   # Sort yarns ###################
   for i,Yarn in enumerate(Yarns):
      M0P=Yarn.GetNode(0).GetPosition()
      M1P=Yarn.GetNode(1).GetPosition()
      V=M1P-M0P
      M2=Yarn.GetNode(2)
      if M2:
         print 'Binder: '+str(i)
         Z.append(Yarn)
      elif CheckX(V):
         print 'Warp: '+str(i)
         X.append(Yarn)
      elif CheckY(V):  
         print 'Weft: '+str(i)
         Y.append(Yarn)
   return X,Y,Z

def GetCentroidsXYZList(YarnList):
   X=[]
   Y=[]
   Z=[]
   for yarn in YarnList:
      Cents=GetCentroidsXYZ(yarn)
      for c in Cents:
         X.append(c.x)
         Y.append(c.y)  
         Z.append(c.z)    
   return X,Y,Z 

def SegIntersect(p1,p2,p3,p4):
   denom = (p4.y-p3.y)*(p2.x-p1.x) - (p4.x-p3.x)*(p2.y-p1.y)
   if not denom:
      print 'Parallel lines!\n'
      return None,None
   dU1 = ((p4.x-p3.x)*(p1.y-p3.y)-(p4.y-p3.y)*(p1.x-p3.x))/denom
   dU2 = ((p2.x-p1.x)*(p1.y-p3.y)-(p2.y-p1.y)*(p1.x-p3.x))/denom
   if dU1<0 or dU1>1 or dU2<0 or dU2>1:
      print 'Unsuccessful!\n'
      return None,None     
   return dU1,dU2     
   
def SortInPlane(YarnSet,Axis,NumInPlane):
   YarnCent=[GetCentroidsXYZ(y) for y in YarnSet]
   YarnHash=[True for y in YarnSet]
   SortedSet=[]
   for i,wa in enumerate(YarnSet[:-1]):
       Layer=[]     
       if YarnHash[i]:
         YarnHash[i]=False 
         Layer.append(YarnCent[i])
         m0=wa.GetNode(0).GetPosition() 
         j=i+1
         for wa1 in YarnSet[i+1:]:
            dirbool=False
            m1=wa1.GetNode(0).GetPosition()
            V=m1-m0             
            #print(V.z)
            if YarnHash[j]:
              #print 'in:'+str(j)+'\n'
              if Axis in ['X','x']:
                 dirbool=TolCheckX(V,0.1)
                 #if dirbool:
                 #print 'Weft yarn'+str(i)+'to '+str(j)+' is:z='+str(V.z)+',x='+str(V.x)+str(dirbool)+'\n'
              elif Axis in ['Y','y']:
                 dirbool=TolCheckY(V,0.1)
                 #if dirbool: 
                 #print 'Warp yarn'+str(i)+'to '+str(j)+' is:z='+str(V.z)+',y='+str(V.y)+str(dirbool)+'\n'
              else:
                 print 'Axis must be either X or Y\n'
                 return 0                          
              if dirbool:
                 Layer.append(YarnCent[j])
                 YarnHash[j]=False  
            else:
               pass     
            j+=1                             
         if Axis in ['X','x']:        
            Layer=sorted(Layer, key=GetX)   
         elif Axis in ['Y','y']:
            Layer=sorted(Layer, key=GetY)
         else:
            print 'Axis must be either X or Y\n'
            return 0              
         SortedSet.append(Layer[:NumInPlane])    
         print([yarn[0].x for yarn in Layer[:NumInPlane]])
         #print([yarn[0].x for yarn in Layer]) 
   return SortedSet

def Find2DIntersection(Yarn1, Yarn2, window):
   mindist=1e5
   #Find closest slave points first 
   for i,slave in enumerate(Yarn1):
       ind=GetClosestPointIndex(Yarn2,slave)
       nearpoint=Yarn2[ind]
       length=GetLength(slave,nearpoint)
       if length<mindist:
          mindist=length
          minind1=i
          minind2=ind
   lim1=len(Yarn1)-1
   lim2=len(Yarn2)-1  
   print(minind1,minind2)    
   X1=[node.x for node in Yarn1]
   Y1=[node.y for node in Yarn1]
   X2=[node.x for node in Yarn2]
   Y2=[node.y for node in Yarn2]
   #Determin neighbouring slave points to consider for local curves using the window variable
   start1=minind1-window
   if start1<0:
      start1=0

   if minind1+window<=lim1:
      end1=minind1+window
   else:
      end1=lim1

   start2=minind2-window
   if start2<0:
      start2=0

   if minind2+window<=lim2:
      end2=minind2+window
   else:
      end2=lim2
   #Fit lines for points neighbouring the closest slave points 
   k1,m1=np.polyfit(np.array(X1[start1:end1]),np.array(Y1[start1:end1]),1)       
   k2,m2=np.polyfit(np.array(X2[start2:end2]),np.array(Y2[start2:end2]),1)

   x10=X1[start1]
   x11=X1[end1]
   dx1=x11-x10
   p11=XY(x10,k1*x10+m1)
   p12=XY(x11,k1*x11+m1)

   x20=X2[start2]
   x21=X2[end2]
   dx2=x21-x20
   p21=XY(x20,k2*x20+m2)
   p22=XY(x21,k2*x21+m2)
   
   dU1,dU2=SegIntersect(p11,p12,p21,p22)
   if dU1 and dU2:
      Point=XY(x10+dU1*dx1,k1*(x10+dU1*dx1)+m1)
   else:
      Point=Yarn1[minind1]
   return Point

def BuildControlMesh(Textile,InPlaneNumWeft,InPlaneNumWarp):
   # Control points are defined a set of points which characterise the architecture pattern distortion
   # The points are supposed to be spaced irregularly because the geometry is distortded
   # By connecting the control point to form a brick mesh of the domain, each element contains
   # a set of points from the point cloud of the measured geometry . 
   # If we assusme that we can "regularise" the grid of points(mesh) by averaging the distances we can
   # come up with an orthogonal brick mesh referring to the corrected geometry 
   # If we assume the same element-wise affine transformations required to get to the regular grid from the 
   # irregular can be applied for the respective point set, a reverse process based on brick element shape functions
   # can be defined. 
   X,Y,Z=YarnTypeSort(Textile)

   window=3
   UpperBinder=[] 
   LowerBinder=[]
   #Replace with the sorting function
   XSort=SortInPlane(X,'Y', InPlaneNumWarp)
   YSort=SortInPlane(Y,'X', InPlaneNumWeft)
   TopWeft=YSort[0]
   UpperWarp=XSort[0]
   NumWeftLay=len(YSort)
   MidWeft=YSort[NumWeftLay//2]
   NumWarpLay=len(XSort)
   UpperMidWarp=XSort[NumWarpLay//2]
   LowerMidWarp=XSort[NumWarpLay//2+1] 
   LowerWarp=XSort[-1]
   BottomWeft=YSort[-1]
   middz=MidWeft[1][0].z   
   
   for yarn in Z:
      if yarn.GetNode(1).GetPosition().z>middz:
         UpperBinder.append(yarn)
      else:
         LowerBinder.append(yarn)   
   upperz=UpperBinder[0].GetNode(1).GetPosition().z + 0.1
   lowerz=LowerBinder[0].GetNode(1).GetPosition().z - 0.1        
   # Top Layer Control Points
   #Control_Points=[[[None for i in range(InPlaneNumWeft-1)] for j in range(InPlaneNumWarp)] for k in range(3)]
   Control_Points_Vec=XYZVector()
   times=0
   for weft in TopWeft:
      for warp in UpperWarp:
         Point1=Find2DIntersection(weft,warp,window)
         #Control_Points[i][j][2]=XYZ(Point1.x,Point1.y,upperz)
         Control_Points_Vec.push_back(XYZ(Point1.x,Point1.y,upperz))
   for weft in MidWeft:
      for j,warp in enumerate(UpperMidWarp):
         Point1=Find2DIntersection(weft,warp,window)
         Point2=Find2DIntersection(weft,LowerMidWarp[j],window)  
         AvPoint=(Point1+Point2)*0.5
         #Control_Points[i][j][1]=XYZ(AvPoint.x,AvPoint.y,middz)
         Control_Points_Vec.push_back(XYZ(AvPoint.x,AvPoint.y,middz))
   for weft in BottomWeft:
      for warp in LowerWarp:
         Point1=Find2DIntersection(weft,warp,window)
         Control_Points_Vec.push_back(XYZ(Point1.x,Point1.y,lowerz))
         #Control_Points[i][j][2]=XYZ(Point1.x,Point1.y,lowerz)
   #print(Control_Points)
   # Find a way to identify cross over points and yarn section in vicinity
   #ControlMesh=XYZVector()

   return Control_Points_Vec

def RegulariseControlPoints(PointVector,PlaneNumX,PlaneNumY):
    NewPointVector=XYZVector()
    D0=XYZ()
    D1=XYZ()
    X=[[p.x,1.0] for p in PointVector]
    Y=[[1.0,p.y] for p in PointVector]
    X=np.array(X)
    Y=np.array(Y)
    kmeansX=KMeans(n_clusters=PlaneNumX,random_state=0).fit(X)
    kmeansY=KMeans(n_clusters=PlaneNumY,random_state=0).fit(Y)
    centX=kmeansX.cluster_centers_[:,0]
    centY=kmeansY.cluster_centers_[:,1] 

    return 0

def PlotScatterProj2D(X,Y,Htitle,Vtitle,Title):
   fig1, ax1=plt.subplots()
   ax1.set_xlabel(Htitle,fontsize=14)
   ax1.set_ylabel(Vtitle,fontsize=14)
   ax1.set_title(Title,fontsize=14)

   ax1.scatter(X,Y, label='Centroid Data')
   #k,m=np.polyfit(np.array(X),np.array(Y),1)
   #ax1.plot(np.array(X),np.array(X)*k+m,color='red', label='Fitted Line')
   #k,m=np.polyfit(np.array(X),np.array(Y),1)
   #ax1.plot(np.array(X),np.array(X)*k+m,color='red')
   ax1.legend(fontsize=14)   
   plt.show()
   return 0

def PlotScatter3D(Vector,ax1,colour):
   #fig1=plt.figure() 
   #ax1=fig1.add_subplot(111,projection='3d')

   X=[p.x for p in Vector]
   Y=[p.y for p in Vector]
   Z=[p.z for p in Vector]
   ax1.scatter(X,Y,Z, marker='.',color=colour)
   ax1.set_xlabel('X',fontsize=14,fontweight='bold')
   ax1.set_ylabel('Y',fontsize=14,fontweight='bold')
   ax1.set_zlabel('Z',fontsize=14,fontweight='bold')
   ax1.set_title('Control Points',fontsize=18,fontweight='bold') 
   #plt.show()
   return 0
def Red(len_ratio):
   if len_ratio>0.5:
      return len_ratio
   else:
      return 0.0
def Blue(len_ratio):
   if len_ratio<0.5:
      return 1.0-len_ratio
   else:
      return 0.0
def Green(len_ratio):
   if len_ratio>=0.5:
      return 2.0-2.0*len_ratio
   else:
      return 2.0*len_ratio                 
def PlotVectorField3D(Vector1,Vector2,ax1):
   #fig1=plt.figure() 
   #ax1=fig1.add_subplot(111,projection='3d')
   
   X1=np.array([p.x for p in Vector1])
   Y1=np.array([p.y for p in Vector1])
   Z1=np.array([p.z for p in Vector1])
   
   X2=np.array([p.x for p in Vector2])
   Y2=np.array([p.y for p in Vector2])
   Z2=np.array([p.z for p in Vector2])
   Lens=np.sqrt(np.power(X2-X1,2)+np.power(Y2-Y1,2)+np.power(Z2-Z1,2))
   max_len=Lens.max()
   min_len=Lens.min()
   NLens=Lens/max_len
   ColourMap=[(Red(l),Green(l),Blue(l)) for l in NLens]
   ax1.quiver(X1,Y1,Z1,X2-X1,Y2-Y1,Z2-Z1,length=0.1, normalize=False, colors=ColourMap)
 
   ax1.set_xlabel('X',fontsize=14)
   ax1.set_ylabel('Y',fontsize=14)
   ax1.set_zlabel('Z',fontsize=14)
   ax1.set_title('Control Points',fontsize=14)    
 
   return min_len,max_len

def PlotProjectedCentroids(Textile):
   Warp, Weft, Binder = YarnTypeSort(Textile)
   #warp
   del Warp[-1] # Incomplete yarns
   del Warp[-1] 
   WarpLeft=[wa for wa in Warp if wa.GetNode(0).GetPosition().y<3.0]
   WarpRight=[wa for wa in Warp if wa.GetNode(0).GetPosition().y>3.0]
   #X,Y,Z=GetCentroidsXYList(Warp)  
   fig1, ax1=plt.subplots()
   ax1.set_xlabel('Y',fontsize=14)
   ax1.set_ylabel('Z',fontsize=14)
   ax1.set_title("Warp Projection",fontsize=14)
   XL,YL,ZL=GetCentroidsXYZList(WarpLeft)
   XR,YR,ZR=GetCentroidsXYZList(WarpRight)
   ax1.scatter(YL+YR,ZL+ZR, label='Centroid Data')
   k,m=np.polyfit(np.array(YL),np.array(ZL),1)
   ax1.plot(np.array(YL),np.array(YL)*k+m,color='red', label='Fitted Line')
   k,m=np.polyfit(np.array(YR),np.array(ZR),1)
   ax1.plot(np.array(YR),np.array(YR)*k+m,color='red')
   ax1.legend(fontsize=14)
   #Weft
   FullWeft=[we for we in Weft if we.GetNode(0).GetPosition().x > 3.5 and we.GetNode(0).GetPosition().x < 7.3 ]
   WeftLeft=[we for we in FullWeft if we.GetNode(0).GetPosition().x<5.0]
   WeftRight=[we for we in FullWeft if we.GetNode(0).GetPosition().x>5.0]   
   XL,YL,ZL=GetCentroidsXYZList(WeftLeft)
   XR,YR,ZR=GetCentroidsXYZList(WeftRight)
   fig2, ax2=plt.subplots()
   ax2.set_xlabel('X',fontsize=14)
   ax2.set_ylabel('Z',fontsize=14)
   ax2.set_title("Weft Projection",fontsize=14)
   ax2.scatter(XL+XR,ZL+ZR, label='Centroid Data')
   k,m=np.polyfit(np.array(XL),np.array(ZL),1)
   ax2.plot(np.array(XL),np.array(XL)*k+m,color='red', label='Fitted Line')
   k,m=np.polyfit(np.array(XR),np.array(ZR),1)
   ax2.plot(np.array(XR),np.array(XR)*k+m,color='red')
   ax2.legend(fontsize=14)
   plt.show()
   
   X,Y,Z=GetCentroidsXYZList(Warp+Weft+Binder)
   PlotScatterProj2D(X,Y,'X','Y','XY Projection')
   return 0

def Extend(Textile,ODomain):
 
    Yarns=Textile.GetYarns()
    Warp, Weft, Binder = YarnTypeSort(Textile)
    D0Y=0
    D0Ycount=0

    OWarp=tuple(Warp)      
    del Warp[-1]
    del Warp[-1]      
    ########################
    ### Compute average translation vectors ################   
    ##  translation vector in Y axis:  
    dysum=0
    Bnum=0
    for i,b in enumerate(Binder[:-1]):
       m0=b.GetNode(2).GetPosition()
       mindist=100.0
       for b1 in Binder[i+1:]:
          m1=b1.GetNode(2).GetPosition()
          V=m1-m0        
          if CheckY(V) and abs(V.y)<=mindist:
             mindist=abs(V.y)
       if mindist<100.0:       
          dysum+=mindist
          Bnum+=1

    Wanum=0
    for i,wa in enumerate(Warp[:-1]):
       m0=wa.GetNode(0).GetPosition()
       mindist=100.0
       for wa1 in Warp[i+1:]:
          m1=wa1.GetNode(0).GetPosition()
          V=m1-m0       
          if V.y>0:
             D0Y+=m0.y## To be updated
             D0Ycount+=1 
          if CheckY(V) and abs(V.y)<=mindist:
             mindist=abs(V.y)
       if mindist<100.0:       
          dysum+=mindist
          Wanum+=1

    dymean=dysum/(Bnum+Wanum)
    print(dymean)    
    D0Ymean=D0Y/D0Ycount

    TVy=XYZ(0.0,2*dymean,0.0)
    #############
    ## Translation vector in X axis ### 
    WeftO=tuple(Weft)

    del Weft[3]#This is a split yarn - not an accurate centroid
    dxsum=0
    Wenum=0
    for i,we in enumerate(Weft[:-1]):
       m0=we.GetNode(0).GetPosition()
       mindist=100.0
       for we1 in Weft[i+1:]:
          m1=we1.GetNode(0).GetPosition()
          V=m1-m0        
          if CheckX(V) and abs(V.x)<=mindist:
             mindist=abs(V.x)
       dxsum+=mindist
       Wenum+=1

    dxmean=dxsum/Wenum
    print(dxmean)    
    
    TVx=XYZ(2*dxmean,0.0,0.0)
    #########################
    NewTextile=CTextile()
    Interpolation=CInterpolationBezier(False, False, False)     
    #### Extend Weft yarns ## 
    tolerance=1.0e-1
    for y in WeftO:    
       
       NewWeft=CYarn()
       NewYarnSection=CYarnSectionInterpPosition()
       backmaster=y.GetNode(0)
       SlaveNodes=list(y.GetSlaveNodes(2))
       SVector=XYZVector()
       SPosList=[]
       SecPtsList=[]
       for s in SlaveNodes:
          SPosList.append(s.GetPosition())
          SecPtsList.append(s.Get2DSectionPoints())
          SVector.push_back(s.GetPosition())
       num=len(SlaveNodes)
       LastSlave=SlaveNodes[-1]
       BackSlavePos=LastSlave.GetPosition()-TVy
       backInd=GetClosestPointIndex(SVector,BackSlavePos)

       if backInd>0:
          for j in range(int(num*0.5)):
            ##Take the next slave node in the list in order to copy the cross-section beyord the current lenght
            
            backnext=SlaveNodes[backInd+j+1]
            #print backnext.GetIndex()
            newfrontpos=backnext.GetPosition()+TVy
            newfrontSecPts=backnext.Get2DSectionPoints()
            SPosList.append(newfrontpos)
            SecPtsList.append(newfrontSecPts)

         # New lead master node
       LastPos=SPosList[-1]
       newlen=(LastPos-backmaster.GetPosition()).y
       for i,pos in enumerate(SPosList):
         Dv=pos-backmaster.GetPosition()
         l=abs(Dv.y/newlen)
         NewYarnSection.AddSection(l,CSectionPolygon(SecPtsList[i]))
       NewWeft.AssignSection(NewYarnSection)
       NewWeft.AddNode(backmaster)
       NewWeft.AddNode(CNode(LastPos))
       NewWeft.AssignInterpolation(Interpolation)
       NewWeft.SetResolution(int(len(SPosList)-1),100) 
       NewTextile.AddYarn(NewWeft)
    for y in OWarp:
       NewTextile.AddYarn(y)
    for y in Binder:
       NewTextile.AddYarn(y)

    D0X=(Weft[1].GetNode(0).GetPosition()).x 
    CP1=XYZ(D0X,D0Ymean,0.0)
    OD0=XYZ()
    OD1=XYZ()
    out=ODomain.GetBoxLimits(OD0,OD1)
    CP2=CP1+TVy+TVx+XYZ(0.0,0.0,OD1.z)
    CDomain=CDomainPlanes(CP1,CP2)
    NewTextile.AssignDomain(CDomain)
    return NewTextile,CDomain

def InterpolateSections(SecPoints0,SecPoints1,coeff):
  InterSecPoints=XYVector()
  for p0 in SecPoints0:
    p1=GetClosestPoint(SecPoints1,p0)
    pI=p0*coeff+p1*(1-coeff)
    InterSecPoints.push_back(pI)
  return InterSecPoints

def EnforcePeriodicity(Textile,Domain):
    Yarns=Textile.GetYarns()
    Warp, Weft, Binder = YarnTypeSort(Textile)
    D0Y=0
    D0Ycount=0
    RepVX=XYZVector()
    RepVY=XYZVector()

    OWarp=tuple(Warp)         
    OD0=XYZ()
    OD1=XYZ()
    out=Domain.GetBoxLimits(OD0,OD1)
    
    NewTextile=CTextile()
    Interpolation=CInterpolationBezier(False, False, False)     
    #### Extend Weft yarns ## 
    tolerance=1.0e-1
    for y in Warp+Weft:
      backmaster=y.GetNode(0)
      frontmaster=y.GetNode(1)
      NewYarn=CYarn()
      NewYarnSection=CYarnSectionInterpPosition()
      SlaveNodes=list(y.GetSlaveNodes(2))
      SVector=XYZVector()
      SPosList=[]
      SecPtsList=[]
      NewSecPtsList=[]
      for s in SlaveNodes:
         SPosList.append(s.GetPosition())
         SecPtsList.append(s.Get2DSectionPoints())
         SVector.push_back(s.GetPosition())
      NewSecPtsList=SecPtsList       
      S0ind=GetClosestPointIndex(SVector,OD0)
      S1ind=GetClosestPointIndex(SVector,OD1)
      DNum=S1ind-S0ind
      Steps=DNum//2
      coeff=0.5
      i=0
      while coeff<=1.0 and i<=Steps:

         NewSecPtsList[S0ind+i]=InterpolateSections(SecPtsList[S0ind+i],SecPtsList[S1ind-i],coeff)     
         NewSecPtsList[S1ind-i]=InterpolateSections(SecPtsList[S1ind-i],SecPtsList[S0ind+i],coeff)
         i+=1
         coeff+=1.0

      LastPos=SPosList[-1]
      newlen=GetLength(LastPos,backmaster.GetPosition())
      for i,pos in enumerate(SPosList):
         Dl=GetLength(pos,backmaster.GetPosition())
         l=abs(Dl/newlen)
         NewYarnSection.AddSection(l,CSectionPolygon(NewSecPtsList[i]))

      NewYarn.AssignSection(NewYarnSection)
      NewYarn.AddNode(backmaster)
      NewYarn.AddNode(frontmaster)
      NewYarn.AssignInterpolation(Interpolation)
      NewYarn.SetResolution(int(len(SPosList)-1)*3,100) 
      NewTextile.AddYarn(NewYarn)   
    for y in Binder:
      NewTextile.AddYarn(y)
    NewTextile.AssignDomain(Domain)
    return NewTextile

def GetPointCloud(Textile):
   #To do 
   Yarns=Textile.GetYarns()
   PointCloud=XYZVector()
   for yarn in Yarns:
      slavenodes=yarn.GetSlaveNodes(2)
      for s in slavenodes:
         secpts=s.GetSectionPoints()
         for point in secpts:
            PointCloud.push_back(point)   
   return PointCloud



def Assign8thDomain(Textile): #under construction
   
   return Textile

class Yarns: # This class, after initialisation, takes sections using InsertSection and either creates a new yarn in the tree or updates an existing one

   def __init__(self,Index):
      self.Index=Index
      self.Type=None
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
         mx=(m0x+m1x)/2
         my=(m0y+m1y)/2
         if self.Sections[S].Direction in ['X','x']:
            self.Type='Warp'
            m0=np.array([MinS.Slice,mx,my])
            m1=np.array([MaxS.Slice,mx,my])
            t=np.array([1.0,0.0,0.0])
            up=np.array([0.0,0.0,1.0])
            M0=MasterNode(0,m0,0.0,t,up)
            M1=MasterNode(1,m1,0.0,t,up)    
         elif self.Sections[S].Direction in ['Y','y']:
            self.Type='Weft'
            m0=np.array([mx,MinS.Slice,my])
            m1=np.array([mx,MaxS.Slice,my])
            t=np.array([0.0,1.0,0.0])
            up=np.array([0.0,0.0,1.0])
            M0=MasterNode(0,m0,0.0,t,up)
            M1=MasterNode(1,m1,0.0,t,up)
         elif self.Sections[S].Direction in ['Z','z']:
            self.Type='Binder'
            if self.Sections[S].Sign:
               if S<0:
                 x=m1x
                 y=m1y
               else:
                 x=m0x
                 y=m0y                    
               m0=np.array([x,y,MinS.Slice])
               m1=np.array([x,y,MaxS.Slice])
               t=np.array([0.0,0.0,1.0])              
               up=np.array([0.0,-1.0,0.0])                 
               M0=MasterNode(0,m0,0.0,t,up)
               M1=MasterNode(1,m1,0.0,t,up) 
            else:
               if S<0:
                 x=m1x
                 y=m1y
               else:
                 x=m0x
                 y=m0y 
                               
               m0=np.array([x,y,DS[2]-MinS.Slice])
               m1=np.array([x,y,DS[2]-MaxS.Slice])              
               t=np.array([0.0,0.0,-1.0])              
               up=np.array([0.0,-1.0,0.0])              
               M1=MasterNode(0,m1,0.0,t,up)
               M0=MasterNode(1,m0,0.0,t,up)    
                                         
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
           #print(NodeList[2*i].Index,NodeList[2*i+1].Index)
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
  cwdvf=cwd+'\\VF55'
  os.chdir(cwdvf)
  
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
  DatFold='\\Data8'
  os.chdir(cwdvf+DatFold)
  # Store Files:
  FileList=[(f.replace('.dat','')).split('_') for f in listdir(cwdvf+DatFold)] # list of info from names
  FileNames=[f for f in listdir(cwdvf+DatFold)] # full names
  # Initialise auxiliary class: 
  MyYarns=Yarns(0)
  i=0
  for file in FileList:
    # File name structure : Direction_YarnIndex_Index_Slice ex. : X_2_1_230 + _-1 depending on the partition sequence. 
    Direction=file[0]
    YarnIndex=int(file[1])
    Index=int(file[2])
    Slice=int(file[3])*ImgRes
    Sign=None
    Polygon=np.genfromtxt(cwdvf+DatFold+'\\'+FileNames[i])*ImgRes
    #Polygon=RotateArray(Polygon,0.5)
    if Direction in ['Z','z']:
       try: 
          Sign=int(file[4])
          Slice=DS[2]-Slice #Reversed Z axis
       except IndexError:
          #Slice=DS[2]-Slice
          pass
       Polygon=Polygon*np.array([1.0,-1.0])+np.array([0.0, DS[1]]) #DS[1]])#([DS[0],DS[1]])
    elif Direction in ['X','x']:
       #Slice=DS[0]-Slice
       Polygon=Polygon*np.array([-1.0,-1.0])+np.array([DS[1],DS[2]])
    elif Direction in ['Y','y']:
       Polygon=Polygon*np.array([1.0,-1.0])+np.array([0.0,DS[2]])#DS[2]])             
    #Populate trees   
    MySection=Section(Index,Slice,Polygon,Direction,Sign)
    MyYarns.InsertSection(YarnIndex,MySection)
    i+=1
  # Add master nodes to join sections and compute final global positions  
  #MyYarns.PrintYarnTree()
  
  MyYarns.AddMasterNodes()
  
  #######
  # TexGen classes initialisation: 
  Textile=CTextile()
  Interpolation=CInterpolationBezier(False, False, False)
  
  #Traverse auxiliary yarn tree and extract yarns
  MyYarnDict=MyYarns.ExtractYarns()
  #Indices used for binder yarns need to be included in a list
  #This helps with identifying and applying the appropriate resolution 
  BinderIndexList=[]#Initialised and populated with yarns with number of partitions>1

  # Iterate yarn class to extact data and populate corresponding TexGen classes
  for y in MyYarnDict:
    MyYarn=MyYarnDict[y]
    MyNodes=MyYarn.Nodes
    NodeList=MyNodes.GetList([])
######  Add extra nodes between partition links 
#    num=len(NodeList)
#    if num>2:
#      for i in range(num//2-1):
#        N0=NodeList[3*i+1]
#         N1=NodeList[3*i+2]
#         dv=N1.Position-N0.Position
#         NmidPos=N0.Position+dv*0.5
#         tan_i=dv/np.linalg.norm(dv)
#         up_i=np.cross(tan_i,np.array([0.0,1.0,0.0]))
#         Nmid=MasterNode(3*i+2,NmidPos,0.0,tan_i,up_i)
#         NodeList.insert(3*i+2,Nmid)
############################# 

    n0y=NodeList[0].Position[1]
    NumSlices=MyYarn.CountSlices(0)
    MySections=MyYarn.Sections
    BinderBool=False
    if len(MySections)>1:
      BinderIndexList.append(y)
      BinderBool=True
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
        
        #Append nodes for translation vector computations:
        #if BinderBool and sec==0:
        #   if 

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
    # Adjusted resolution depending on yarn type. 
    if BinderBool:
      CYarn0.SetResolution(int(NumSlices*0.6),50)
    else:
      CYarn0.SetResolution(int(NumSlices*0.8),100) 
    Textile.AddYarn(CYarn0)

  #Save tg3 file  
  
  #NewTex,NewDomain=Extend(Textile,CDomain)
  #NewTex=EnforcePeriodicity(NewTex,NewDomain)
  #NewTex.AssignDomain(CDomain)
  ControlPts=BuildControlMesh(Textile,CDomain,4,2)
  DataPts=GetPointCloud(Textile)
  mesh,regmesh=cm.BuildMesh(ControlPts,4,2,3)

  #os.chdir(cwd)
  #file=open('TestMesh.inp','w')
  #file.write('*Node\n')
  #for i in gnodes:
  #    node=gnodes[i]
  #    pos=node.Position
  #    file.write(str(i)+', '+str(pos.x)+', '+str(pos.y)+', '+str(pos.z)+'\n')
  #file.write('*Element, Type=C3D8R\n')
  #for el in mesh:
  #   string=''
  #   for n in el.Connectivity:
  #      string+=', '+str(n)
  #   file.write(str(el.Index)+string+'\n')
  #file.write('*NSet, NSet=hold, Unsorted\n1\n*Boundary\nhold,1,1\n')
  #file.write('*Step, Name=myStep\n*Static\n*End Step\n')
#
  #file.close()   
  #
  
  regmesh=cm.UndistortedPointCloud(DataPts,mesh,regmesh)
  fig=plt.figure()
  font={'family':'normal',
        'weight':'bold',
        'size'  : 14}
  plt.rc('font',**font)      
  ax=fig.add_subplot(111,projection='3d')
  PlotScatter3D(ControlPts,ax,'black')
  
  #Plot Region : Comment lines according to designated plotting sections

  ###Plot scatter of data points in control mesh
  #numofel=len(mesh)
  #CoMap=[(Red(r),Green(r),Blue(r)) for r in np.linspace(0.0,1.0,numofel)]
  #for i,e in enumerate(mesh):
  #   PlotScatter3D(e.DataPoints,ax,CoMap[i])
  #ax.set_xlim(0.0,10.0)
  #ax.set_ylim(0.0,10.0)
  #ax.set_zlim(0.0,10.0)
  #ax.set_aspect('auto',adjustable=None)
  #ax.set_axis_off()
  #norm=mpl.colors.Normalize(vmin=0,vmax=int(numofel))
  #cmap=mpl.colors.ListedColormap(CoMap,name='MyMap')     
  #SM=mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
  #SM.set_array([])
  #fig.colorbar(SM, label='Element')  
  ### end
  
  ###Plot displacement vector field:
  #InPoints=XYZVector()
  #OutPoints=XYZVector()
  #for i in range(len(regmesh)):
  #    for j in range(len(regmesh[i].DataPoints)):
  #       InPoints.push_back(mesh[i].DataPoints[j])
  #       OutPoints.push_back(regmesh[i].DataPoints[j]) 
  #minlen,maxlen=PlotVectorField3D(InPoints,OutPoints,ax)
  ax.set_xlim(0.0,9.0)
  ax.set_ylim(0.0,9.0)
  ax.set_zlim(0.0,6.0)
  ax.set_aspect('auto',adjustable=None)
  #CoMap=[(Red(r),Green(r),Blue(r)) for r in np.linspace(0.0,1.0,20)]
  #norm=mpl.colors.Normalize(vmin=minlen,vmax=maxlen)
  #cmap=mpl.colors.ListedColormap(CoMap,name='MyMap')
  #SM=mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
  #SM.set_array([])
  #fig.colorbar(SM, label='Displacement [mm]')
  ###end
  
  plt.show() # Comment accordingly

  #Finalise and save TG3 file:
  #Textile.AssignDomain(CDomain)
  #AddTextile('Rec9',NewTex)
  #AddTextile('Rec9',Textile)
  #SaveToXML(cwd+'\\Reconstruction10.tg3',"",OUTPUT_STANDARD)


    




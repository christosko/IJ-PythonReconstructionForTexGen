from cmath import isnan
from turtle import position
import numpy as np
import sys#,os
from sklearn.cluster import KMeans
#workdir=os.getcwd()
#Requires a minimum build of texgen and a path pointing to Python libraries
TexGenPath='C:\\Python27\\Lib\\site-packages\\TexGen'
sys.path.append(TexGenPath)
from TexGen.Core import*

class BrickElement:
    def __init__(self,Index):
        self.Index=int(Index)
        self.Nodes=[]
        self.Connectivity=[]

    def InsertNode(self,node):
        if isinstance(node,Node):
           self.Nodes.append(node)
           return 0
        else:
           print 'Invalide type for node\n'
           return None

class Node:
     def __init__(self,Index,Position):
        self.Index=int(Index)
        if isinstance(Position,XYZ):
           self.Position=Position
        else:
           self.Position=None  
           print 'Invalid position type - Needs to be TexGen.Core.XYZ()\n'             
        self.ElementsShared=[]
     def Move(self,Vector):
        if isinstance(Vector,XYZ):
           self.Position+=Vector
           return True
        else:
           print 'Invalid translation vector type: Must be TexGen.Core.XYZ()\n'
           return False

class GlobalShapeFunctions:
    def __init__(self,Element):
        if isinstance(Element,BrickElement):
           Nodes=Element.Nodes
           if Nodes:           
              p=[]
              ShapeCoeffs=np.empty[8,8]
              G=np.empty([8,8])
              for i,n in enumerate(Nodes):
                 p=n.Position
                 Grow=np.array([1,p.x,p.y,p.z,p.x*p.y,p.x*p.z,p.y*p.z,p.x*p.y*p.z])
                 G[i]=Grow
              InvG=np.linalg.inv(G) 
              for i,n in enumerate(Nodes):
                 delta=np.zeros([1,8])
                 delta[0][i]=1
                 coeffs=np.dot(InvG,delta)
                 ShapeCoeffs[i]=coeffs
              self.Coefficients=ShapeCoeffs 
           else:
              print 'Failed to initialise: Nodes not found\n' 
              self.Coefficients=np.zeros([8,8])            
        else:
           print 'Failed to intialise: Void element\n'
           self.Coefficients=np.zeros([8,8])  

    def Value(self,i,position):
        if i in range(1,8):
           if isinstance(position,XYZ):
              x=position.x
              y=position.y
              z=position.z
              Grow=np.array([1,x,y,z,x*y,x*z,y*z,x*y*z])
              ci=self.Coefficients[i]
              ciT=np.transpose(ci)
              return np.dot(Grow,ciT)
           else:
              print 'Invalid position type - Needs to be TexGen.Core.XYZ()'
              return None   
        else:
           print 'i must be an integer from 0 to 7'
           return None    

class LocalShapeFunctions:
    def __init__(self):
        self.LocalNodeCo=np.array([[-1,-1,-1],
                                   [1,-1,-1],  
                                   [1,1,-1],
                                   [-1,1,-1],
                                   [-1,-1,1],
                                   [1,-1,1],  
                                   [1,1,1],
                                   [-1,1,1]])
    def Value(self,i,position):
        if i in range(1,8):
           if isinstance(position,np.array):
              zeta=position[0]
              eta=position[1]
              ksi=position[2]
              nodeCo=self.LocalNodeCo[i]
              zi=nodeCo[0]
              hi=nodeCo[1]
              fi=nodeCo[2]
              eighth=1/8
              return eighth*(1+zi*zeta)*(1+hi*eta)*(1+fi*ksi)
           else:
              print 'Invalid position type - Needs to be TexGen.Core.XYZ()'
              return None   
        else:
           print 'i must be an integer from 0 to 7'
           return None    

def SortPlaneLabels(Labels,CCentres):
    sortindices=CCentres[:,0].argsort()
    labelswap={}
    for i,ind in enumerate(sortindices):
        labelswap[ind]=i
    SortedLabels=[labelswap[i] for i in Labels]    
    return SortedLabels

def BuildMesh(Points,NumYZPlanes,NumXZPlanes,NumXYPlanes):

    X=[[p.x,0.0] for p in Points]
    Y=[[p.y,0.0] for p in Points]
    Z=[[p.z,0.0] for p in Points]
    X=np.array(X)
    Y=np.array(Y)
    Z=np.array(Z)
    kmeansX=KMeans(n_clusters=NumYZPlanes,random_state=0).fit(X)
    kmeansY=KMeans(n_clusters=NumXZPlanes,random_state=0).fit(Y)
    kmeansZ=KMeans(n_clusters=NumXYPlanes,random_state=0).fit(Z)
    xmax=NumYZPlanes
    ymax=NumXZPlanes
    zmax=NumXYPlanes
    xlabels=SortPlaneLabels(kmeansX.labels_,kmeansX.cluster_centers_)
    ylabels=SortPlaneLabels(kmeansY.labels_,kmeansY.cluster_centers_)
    zlabels=SortPlaneLabels(kmeansZ.labels_,kmeansZ.cluster_centers_)
    
    GlobalNodes={}
    
    for i,x in enumerate(xlabels):
      index=x+1+xmax*ylabels[i]+xmax*ymax*zlabels[i]  
      GlobalNodes[index]=Node(index,Points[i])

    addY=xmax
    addZ=xmax*ymax
    Conn=[[i+j*addY+k*addZ, 
      1+i+j*addY+k*addZ,
      1+i+(j+1)*addY+k*addZ,      
      i+(j+1)*addY+k*addZ,      
      i+j*addY+(k+1)*addZ, 
      1+i+j*addY+(k+1)*addZ, 
      1+i+(j+1)*addY+(k+1)*addZ,      
      i+(j+1)*addY+(k+1)*addZ,         
      ] for i in range(1,xmax) for j in range(ymax-1) for k in range(zmax-1)]

    Elements=[]
    for i,e in enumerate(Conn):
       myElement=BrickElement(i+1)
       for n in e:
         myElement.InsertNode(GlobalNodes[n])
         GlobalNodes[n].ElementsShared.append(i)
       myElement.Connectivity=e  
       Elements.append(myElement)        
    return Elements,GlobalNodes

if __name__=='__main__':
    test=XYZVector()
    for i in range(5):
        for j in range(4):
            for k in range(3):
                test.push_back(XYZ(i*2.6,j*1.8,k*1.3))
    mesh=BuildMesh(test,5,4,3)
    file=open('TestMesh.inp','w')
    file.write('*Node\n')
    for i,node in enumerate(test):
        file.write(str(i)+', '+str(node.x))



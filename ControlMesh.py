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
              ShapeCoeffs=np.empty([8,8])
              G=np.empty([8,8])
              for i,n in enumerate(Nodes):
                 p=n.Position
                 Grow=np.array([1,p.x,p.y,p.z,p.x*p.y,p.x*p.z,p.y*p.z,p.x*p.y*p.z])
                 G[i]=Grow
              InvG=np.linalg.inv(G) 
              for i,n in enumerate(Nodes):
                 delta=np.zeros([1,8])
                 delta[0][i]=1.0
                 coeffs=np.dot(InvG,np.transpose(delta))
                 ShapeCoeffs[i]=np.transpose(coeffs)
              self.Coefficients=ShapeCoeffs 
           else:
              print 'Failed to initialise: Nodes not found\n' 
              self.Coefficients=np.zeros([8,8])            
        else:
           print 'Failed to intialise: Void element\n'
           self.Coefficients=np.zeros([8,8])  

    def Value(self,i,position):
        if i in range(8):
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
        self.LocalNodeCo=np.array([[-1.0,-1.0,-1.0],
                                   [1.0,-1.0,-1.0],  
                                   [1.0,1.0,-1.0],
                                   [-1.0,1.0,-1.0],
                                   [-1.0,-1.0,1.0],
                                   [1.0,-1.0,1.0],  
                                   [1.0,1.0,1.0],
                                   [-1.0,1.0,1.0]])
    def Value(self,i,position):
        if i in range(8):           
           zeta=position[0]
           eta=position[1]
           ksi=position[2]
           nodeCo=self.LocalNodeCo[i]
           zi=nodeCo[0]
           hi=nodeCo[1]
           fi=nodeCo[2]
           eighth=1.0/8.0
           return eighth*(1.0+zi*zeta)*(1.0+hi*eta)*(1.0+fi*ksi) 
        else:
           print 'i must be an integer from 0 to 7'
           return None    

class BrickElement:
    def __init__(self,Index):
        self.Index=int(Index)
        self.Nodes=[]
        self.NodePositions=XYZVector()
        self.Connectivity=[]
        self.ShapeFunctions=[]
        self.DataPoints=XYZVector()
        self.IndexMatch={}

    def InsertNode(self,node):
        if isinstance(node,Node):
           self.Nodes.append(node)
           self.NodePositions.push_back(node.Position)
           return 0
        else:
           print 'Invalide type for node\n'
           return None

    def InsertPoint(self,Point,GlobalIndex=0):
        if isinstance(Point,XYZ):
           ind=len(self.DataPoints)
           self.IndexMatch[ind]=GlobalIndex 
           self.DataPoints.push_back(Point)
           return 0
        else:
           pass 
           print 'Invalid point type - Needs to be TexGen.Core.XYZ()\n'
           return None

    def InitGlobalShapeFunctions(self):
        self.ShapeFunctions=[GlobalShapeFunctions(self) for i in range(8)]
        return 0        

    #def GetGlobalShapeValue(self,node_index,point_position):
    #    return self.ShapeFunctions[node_index].Value(node_index,point_position)


def SortPlaneLabels(Labels,CCentres):
    #Labels of coplanar cluster centres are off - this sorts them
    sortindices=CCentres[:,0].argsort()
    labelswap={}
    for i,ind in enumerate(sortindices):
        labelswap[ind]=i
    SortedLabels=[labelswap[i] for i in Labels]    
    return SortedLabels,np.array([CCentres[i] for i in sortindices])

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
    xclust=kmeansX.cluster_centers_
    yclust=kmeansY.cluster_centers_
    zclust=kmeansZ.cluster_centers_
    xlabels,xclust=SortPlaneLabels(kmeansX.labels_,xclust)
    ylabels,yclust=SortPlaneLabels(kmeansY.labels_,yclust)
    zlabels,zclust=SortPlaneLabels(kmeansZ.labels_,zclust)

    GlobalNodes={}
    RegGlobalNodes={}
    for i,x in enumerate(xlabels):
      index=x+1+xmax*ylabels[i]+xmax*ymax*zlabels[i]  
      GlobalNodes[index]=Node(index,Points[i])
      regpos=XYZ(xclust[xlabels[i]][0],yclust[ylabels[i]][0],zclust[zlabels[i]][0])
      RegGlobalNodes[index]=Node(index,regpos)

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
    RegElements=[]
    for i,e in enumerate(Conn):
       myElement=BrickElement(i+1)
       myRegElement=BrickElement(i+1)
       for n in e:
         myElement.InsertNode(GlobalNodes[n])
         myRegElement.InsertNode(RegGlobalNodes[n])
         GlobalNodes[n].ElementsShared.append(i)
         RegGlobalNodes[n].ElementsShared.append(i)
        
       myElement.Connectivity=e  
       myRegElement.Connectivity=e
       myElement.InitGlobalShapeFunctions()
       
       Elements.append(myElement) 
       RegElements.append(myRegElement)       
    return Elements, RegElements, GlobalNodes

def TrasformationMatrix(FaceNodes):
    #A matrix used to epxress data points in element's face local coordinates to identify inclusion
    #Nodes are assumed to be coplanar. Needs checking 
    #Face nodes are in counter clockwise order (right hand rule)
    P1=FaceNodes[0]
    P2=FaceNodes[1]
    P3=FaceNodes[2]
    P4=FaceNodes[3]
    v12=P2-P1
    v14=P4-P1
    outv=CrossProduct(v12,v14)
    w=outv*(1/GetLength(outv))
    #print(w)
    u=v12*(1/GetLength(v12))
    v=CrossProduct(w,u)
    i=XYZ(1.0,0.0,0.0)
    j=XYZ(0.0,1.0,0.0)
    k=XYZ(0.0,0.0,1.0)
    T=np.array([[DotProduct(u,i),DotProduct(u,j),DotProduct(u,k)],
    [DotProduct(v,i),DotProduct(v,j),DotProduct(v,k)],
    [DotProduct(w,i),DotProduct(w,j),DotProduct(w,k)]])

    return T 

def PointInHex(Point,NodesVector):
    # Node list is in abaqus standard order : bottom 
    # face counter clockwise and top face counter clockwise
    # 6 faces with outward pointing normal vectors must be defined    
    IN=[False for i in range(6)]
    N=NodesVector
    #Faces with node orders subject to right hand rule for outward pointing normal vector
    FaceNodesInds=[
    [0,3,2,1],
    [0,4,7,3],
    [4,5,6,7],
    [5,1,2,6],
    [0,1,5,4],
    [3,7,6,2]    
    ]
    FaceNodes=[[N[i] for i in n] for n in FaceNodesInds]

    TMats=[TrasformationMatrix(Face) for Face in FaceNodes]

    PointV=np.array([Point.x,Point.y,Point.z])
    for i,T in enumerate(TMats):
        PointVT=np.dot(T,PointV)
        Node1Pos=np.array([FaceNodes[i][0].x,FaceNodes[i][0].y,FaceNodes[i][0].z])
        N1T=np.dot(T,Node1Pos)
        LPointVT=PointVT-N1T
        if LPointVT[2]<=0:
            IN[i]=True
  
    return all(IN)

def UndistortedPointCloud(PointCloud,DistMesh,RegMesh):
    
    # Points need to be labeled based on the element the lie in 
    # This is achieved by checking if a point lies in the hexahedron 
    Hash=[True for point in PointCloud]
    #MemUPC=[None for point in PointCloud]
    UPointCloud=XYZVector(len(PointCloud))
    for i,p in enumerate(PointCloud):
        UPointCloud[i]=p 
    LSF=[LocalShapeFunctions() for i in range(8)]
    nodezeta=[-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0]
    nodeeta=[-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0]
    nodeksi=[-1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0]    
    for ind,element in enumerate(DistMesh):
        nodeset=element.NodePositions
        StoreInds=[]
        distpoints=XYZVector()
        for j,p in enumerate(PointCloud):
            if Hash[j]:
                if PointInHex(p,nodeset):
                   element.InsertPoint(p,j)
                   distpoints.push_back(p)
                   StoreInds.append(j)
                   Hash[j]=False

        #distpoints=element.DataPoints
        regnodepos=RegMesh[ind].NodePositions
        SF=element.ShapeFunctions

        for pi,p in enumerate(distpoints):
            pz=0.0
            pe=0.0
            pk=0.0
            for i in range(8):
                pz+=nodezeta[i]*SF[i].Value(i,p)
                pe+=nodeeta[i]*SF[i].Value(i,p)
                pk+=nodeksi[i]*SF[i].Value(i,p)
            newx=0.0
            newy=0.0
            newz=0.0
            for i in range(8):
                newx+=regnodepos[i].x*LSF[i].Value(i,np.array([pz,pe,pk]))
                newy+=regnodepos[i].y*LSF[i].Value(i,np.array([pz,pe,pk]))
                newz+=regnodepos[i].z*LSF[i].Value(i,np.array([pz,pe,pk]))    
             
            UPointCloud[StoreInds[pi]]=XYZ(newx,newy,newz)
            #print(UPointCloud[StoreInds[pi]]-p)   
            RegMesh[ind].InsertPoint(XYZ(newx,newy,newz))
            
        #RegMesh[ind].IndexMatch=element.IndexMatch
    #UPointCloud=XYZVector()
    #for el in RegMesh:
    #
    #    for i,p in enumerate(el.DataPoints):
    #        MemUPC[el.IndexMatch[i]]=p
            
    #for i,p in enumerate(MemUPC):
    #    try: 
    #      UPointCloud.push_back(p)
    #    except ValueError:
    #      print('None value in position: '+str(i))    
    
    return RegMesh, UPointCloud

if __name__=='__main__':
    #test=XYZVector()
    #for i in range(5):
    #    for j in range(4):
    #        for k in range(3):
    #            test.push_back(XYZ(i*2.6,j*1.8,k*1.3))
    #mesh,gnodes=BuildMesh(test,5,4,3)
    #file=open('TestMesh.inp','w')
    #file.write('*Node\n')
    #for i,node in enumerate(test):
    #    file.write(str(i)+', '+str(node.x))
    N=XYZVector()
    N.push_back(XYZ(0,0,0))
    N.push_back(XYZ(1,0,0))#(1.2,-0.4,-0.3))
    N.push_back(XYZ(1,1,0))#(0.8,1.3,0.0))
    N.push_back(XYZ(0,1,0))#(-0.1,1.34,-0.1))
    N.push_back(XYZ(0,0,1))#(0.08,-0.1,1.04))
    N.push_back(XYZ(1,0,1))#(1.1,0.05,0.9))
    N.push_back(XYZ(1,1,1))#(1.23,1.13,1.03))
    N.push_back(XYZ(0,1,1))#(-0.04,0.94,1.1))
    P=XYZ(0.999,0.5,0.5)
    print(PointInHex(P,N))




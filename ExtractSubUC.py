def Extend(Textile,ODomain):
 
    #Yarns=Textile.GetYarns()
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
    TVx=XYZ(2*dxmean,0.0,0.0)
    #########################
    NewTextile=CTextile()
    Interpolation=CInterpolationBezier(False, False, False)     
    #### Extend Weft yarns ## 

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
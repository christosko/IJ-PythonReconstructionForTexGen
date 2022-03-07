import os
import sys
from os import listdir
import numpy as np
#Get in the appropriate VF folder
cwd=os.getcwd() 
ModelLocation=cwd+'\\VF64'
os.chdir(ModelLocation)

DataLocation=ModelLocation+'\\Data2'
FileNames=[f for f in listdir(DataLocation)] # full names
FileList=[(f.replace('.dat','')).split('_') for f in listdir(DataLocation)]
for i,file in enumerate(FileList):   
   Direction=file[0]
   YarnIndex=file[1]
   SectionIndex=file[2]
   Slice=file[3]
   Polygon=np.genfromtxt(DataLocation+'\\'+FileNames[i])
   OldFile=open(DataLocation+'\\'+FileNames[i],'r')
   Lines=OldFile.readlines()
   OldFile.close()
   
   m0x=np.sum(Polygon[:,0])/len(Polygon[:,0])
   m0y=np.sum(Polygon[:,1])/len(Polygon[:,1])
   
   if YarnIndex=='5' and m0y<300:
       print(YarnIndex,m0y)
       Newfile=open(DataLocation+'\\'+Direction+'_'+'6_'+SectionIndex+'_'+Slice+'.dat','w')
       Newfile.writelines(Lines)
       Newfile.close()
       os.remove(DataLocation+'\\'+FileNames[i]) 
   elif YarnIndex=='6' and m0y>300:
       print(YarnIndex,m0y)
       Newfile=open(DataLocation+'\\'+Direction+'_'+'5_'+SectionIndex+'_'+Slice+'.dat','w')
       Newfile.writelines(Lines)
       Newfile.close()     
       os.remove(DataLocation+'\\'+FileNames[i])       
   
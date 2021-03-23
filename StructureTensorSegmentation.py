import numpy as np 
import os 
import sys

cwd=os.getcwd()


class Image:
	def __init__(self,Array):
		self.I=Array
		self.Voxels=VoxelTree

	def StructureTensor(self,SubDomain,Window):
        D=np.array([1,-8,0,8,-1])
        for i in range(len(self.I[:,1,1])):
        	for j in range(len(self.I[1,:,1])):
        		for k in range(len(self.I[1,1,:])):
        			dIdx=D.dot(I[i-2:i+2,j,k])/12
        			dIdy=D.dot(I[i,j-2:j+2,k])/12
        			dIdz=D.dot(I[i,j,k-2:k+2])/12
        			s=np.array([[pow(dIdx,2),dIdx*dIdy,dIdx*dIdz],
        				        [dIdy*dIdx, pow(dIdy,2), dIdy*dIdz],
        				        [dIdz*dIdx, dIdz*dIdy, pow(dIdz,2)]])
                    S=np.empty(3,3)
                    for x in range(i-2,i+2):
                    	for y in range(j-2,j+2):
                    		for z in range(k-2,k+2):



 

	    return ST	
    def PartialDerivative(self,I,Dir,h=1,w):
    	for i in range(w,len(I[:,1,1])-w)
        return PD	
    def Eigen(self,ST)   
        return eValues,eVectors     
class Voxel:
	def __init__(self,Pos,G,AG,s,S,)

def kMeans():

ImArray=np.genfromtxt(Name)

Ι=Image(ImArray)

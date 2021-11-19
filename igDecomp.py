import os
import numpy as np 
from numpy import linalg as la
import math as m
ST=np.array([[69.59885116598096, 47.04295267489709],[47.04295267489709, 31.797067901234453]])
w,v=la.eig(ST)

p1=ST[0,0]+ST[1,1]
p2=(ST[0,0]*ST[1,1])-(ST[0,1]**2)
p3=m.sqrt(p1**2-4*p2**2)
l0=p1-p3
#print(p1,p2,p3)

v1,v2,v3=-0.9929523384194845, -0.05843073875569266, -0.10310917706849099

sq=m.sqrt(v1**2+v2**2+v3**2)
print(w)
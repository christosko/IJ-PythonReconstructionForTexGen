import numpy as np

import math as m
import os
import sys
from os import listdir

os.chdir('D:/Code_Repos')
os.system('abaqus python ODB_results.py -o '+'ID00000SM.odb'+' -n ConstraintsDriverXX')
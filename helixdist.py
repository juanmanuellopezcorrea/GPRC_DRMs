#helixsdist

from genericpath import exists
import mdtraj as mdt
import numpy as np
import glob as glob
import pandas as pd
import os as os
import sys
from math import sqrt
import pickle
import pathlib
from math import sqrt
import funtions as fn
# DISTANCES
root="/Users/Juanlopez/Documents/projects/Unsupervised_Methods_Molecular_Dynamics/UnsupervisedMD/dataset"

proteinlst = ["ID115"] # esto esta realizado para un tipo de proteina

lsthelixs=[[32,46,60],[70,83,96],[106,121,136],[150,160,171],[200,215,229],[284,298,312],[322,331,341]]

#
#lsthelixs=[[32,46,60],[70,83,96]]#,[106,121,136],[150,160,171],[200,215,229],[284,298,312],[322,331,341]]






#b=fn.resspoints(root,proteinlst)
a=fn.helix(root,proteinlst,lsthelixs)






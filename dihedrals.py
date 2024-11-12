
#Coordinates

import funtions as fn
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
import glob as glob
import os as os

from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import umap
import glob
import os
import numpy as np
import funtions as fn
# DISTANCES

#root="/Users/Juanlopez/Documents/projects/Unsupervised_Methods_Molecular_Dynamics/UnsupervisedMD/dataset"
 # esto esta realizado para un tipo de proteina

#root="/home/usuaris/jmlopez/Unsup_dataset/" # root cluste



#fn.dihedrals_angles(root,proteinpath)
prot = str("ID114")
datapath="./dataset/"+prot+"/coordinates/angles/*.csv"


UMAP_angles_residues(datapath,prot)



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
# Uniform Manifold Approximation and Projection
#https://umap-learn.readthedocs.io/en/latest/basic_usage.html
#https://github.com/lmcinnes/umap




#datapath="./dataset/"+prot+"/coordinates/an/*.csv"
#fn.extract_angles(datapath,prot)

#fn.UMAP_angles(datapath,prot)


prot = str("ID115")
datapath="./dataset/"+prot+"/coordinates/angles/*.csv"


#fn.UMAP_angles_residues(datapath,prot)
#fn.UMAP_angles(datapath,prot)

#fn.PCA_angles_residues(datapath,prot)
fn.PCA_angles(datapath,prot)





#fn.UMAP_helix(datapath)

#lstrprot=["ID114","ID115"]
#fn.UMAP_3Dposition(lstrprot)
#fn.PCA_3Dposition(lstrprot)
#fn.UMAP_3Dpositionresidues(lstrprot)
#fn.PCA_3Dpositionresidues(lstrprot)

#fn.UMAP_3Dposition(lstrprot)








#datapath="dataset/*/coordinates/*/embeddings/*.csv"
#fn.graph(datapath)

















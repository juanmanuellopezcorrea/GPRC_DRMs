import numpy as np
import pandas as pd
import os
import glob as glob
import funtions as fn
import seaborn as sns

prot=str("ID115")

#datapath="./dataset/"+prot+"/coordinates/angles/*.csv"
#b= PCA_angles(datapath,prot)




#datapath="./dataset/"+prot+"/coordinates/helix/*.csv"
#fn.PCA_helix(datapath)






lstrprot=["ID115","ID114"]
fn.PCA_3Dposition(lstrprot)


datapath="dataset/*/coordinates/*/embeddings/*.csv"
fn.graph(datapath)



datapath="dataset/*/coordinates/*/embeddings/*.csv"


fn.graph(datapath)







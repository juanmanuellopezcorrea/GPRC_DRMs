#Coordinates
import funtions as fn
#scp -r    jmlopez@cluster-rdlab.cs.upc.edu:/home/usuaris/jmlopez/Unsup_dataset/ID114/coordinates/*.csv  /Users/Juanlopez/Documents/projects/Unsupervised_Methods_Molecular_Dynamics/UnsupervisedMD/dataset/ID114/coordinates

# DISTANCES

root="/Users/Juanlopez/Documents/projects/Unsupervised_Methods_Molecular_Dynamics/UnsupervisedMD/dataset"
proteinpath= ["ID114","ID115"]

fn.resspoints(root,proteinpath)


#lstrprot= ["ID114","ID115"]
#fn.UMAP_3Dposition(lstrprot)
#fn.PCA_3Dposition(lstrprot)

#angles
#fn.dihedrals_angles(root,proteinpath)

import funtions as fn
prot = str("ID115")
datapath="./dataset/"+prot+"/coordinates/angles/*.csv"
#fn.extract_angles(datapath,prot)
fn.UMAP_angles(datapath,prot)
#fn.PCA_angles(datapath,prot)

prot = str("ID115")
#datapath="./dataset/"+prot+"/coordinates/angles/*.csv"
#fn.extract_angles(datapath,prot)
fn.UMAP_angles(datapath,prot)
fn.PCA_angles(datapath,prot)

#helix
import funtions as fn
root="/Users/Juanlopez/Documents/projects/Unsupervised_Methods_Molecular_Dynamics/UnsupervisedMD/dataset"
proteinlst = ["ID114","ID115"] #
lsthelixs=[[32,46,60],[70,83,96],[106,121,136],[150,160,171],[200,215,229],[284,298,312],[322,331,341]]
#fn.helix(root,proteinlst,lsthelixs)




import funtions as fn
prot = str("ID115")
datapath="./dataset/"+prot+"/coordinates/helix/*.csv"
#fn.UMAP_helix(datapath)
#fn.PCA_helix(datapath)
fn.TICA_helix(datapath)

import funtions as fn
prot = str("ID114")
datapath="./dataset/"+prot+"/coordinates/helix/*.csv"
#fn.helix (root, proteinlst, lsthelixs)
#fn.UMAP_helix(datapath)
#fn.PCA_helix(datapath)
fn.TICA_helix(datapath)

#scartter
import funtions as fn
datapath="./dataset/*/coordinates/3Drespositions/embeddings/*/*/*.csv"
fn.graph(datapath)

#ver estos
datapath="./dataset/*/coordinates/angles/embeddings/*/*/*.csv"
fn.graph(datapath)


import funtions as fn
datapath="./dataset/*/coordinates/helix/embeddings/*.csv"
fn.graph(datapath)


#Coordinates
from genericpath import exists
import mdtraj as mdt
import numpy as np
import glob as glob
import pandas as pd
import os as os
import glob as glob
import sys
from math import sqrt
import pickle
import pathlib
# DISTANCES

root="/Users/Juanlopez/Documents/projects/Unsupervised_Methods_Molecular_Dynamics/UnsupervisedMD/dataset"
#proteinpath = "ID114"  # esto esta realizado para un tipo de proteina
#nsteps=30
#root="/home/usuaris/jmlopez/Unsup_dataset/" # root cluster


#variable parameters
#proteinpath = str(sys.argv[1])# "ID114" ID115  # esto esta realizado para un tipo de proteina
#nsteps=int(sys.argv[2])

proteinpath= "ID114"

nsteps=50

tjs=glob.glob(os.path.join(root,proteinpath,"*.xtc"))

pdb = glob.glob(os.path.join(root,proteinpath,"*.pdb"))
for xtc in tjs[:]:#TODO quitar
    nframesin=0
    idframe=2490
    dicframes = {"idframe":[], "idatom":[], "massatom":[],"idres":[],"nameres":[], "x":[],"y":[],"z":[]}
    traj = mdt.load(xtc, top=pdb[0])
    traj=traj.remove_solvent()
    #idframe=0#corregir esto
    for frame in traj[2490:2500]:
        t = frame.topology
        with open(pdb[0]) as f:
            lines = f.readlines()
            #print(lines)
        for i in t.atoms:
            #print(lines[i])
            dicframes["idres"].append(lines[i.index + 4][22:27])
            dicframes["idatom"].append(i.index)
            #dicframes["idres"].append(int (i.residue.index))
            dicframes["massatom"].append(int (i.element.mass))
            dicframes["nameres"].append(lines[i.index + 4][17:22])
            dicframes["idframe"].append(idframe)
        idframe=idframe+1
        for l in frame.xyz:#TODO
                for at in l:
                    dicframes["x"].append(float(at[0]))
                    dicframes["y"].append(float(at[1]))
                    dicframes["z"].append(float(at[2]))
    df_features=pd.DataFrame(dicframes)

    #alfa carbons
    df_features=df_features.loc[df_features.nameres.astype(str)==str("CA")] # only alfa carbons

d
    nframesin=nframesin+nsteps
    #Residues
    ideresidues = df_features.idres.unique()  #ver de eliminar los POPC
    frs = df_features.idframe.unique().tolist()
    dicres={"idress":[],"massx":[],"massy":[],"massz":[],"idframe":[]}

    for i in ideresidues[:]:
        # residues
        df_1 = df_features.loc[df_features.idres.astype(str) == str(i)]
        for f in frs:
            #print(f)
            df_n = df_1.loc[(df_1.idframe.astype(str) == str(f))]
            # centro de masa
            #dfdm = (df_n[["x", 'y', 'z']].multiply(df_n["massatom"], axis="index").sum()) / (
                #df_n.massatom.sum())  # mass center of fn

            dicres["massx"].append(dfdm.x)
            dicres["massy"].append(dfdm.y)
            dicres["massz"].append(dfdm.z)
            dicres["idframe"].append(f)
            dicres["idress"].append(i)
            #dicres["ressname"].append(df_1.nameres.unique()[0])
    df = pd.DataFrame(dicres)


    foldercoord=os.path.join(root,proteinpath,"coordinates" )

    name=xtc.split("/")[-1].split(".xtc")[0]+"_"+str("2500")+"_"+"res.csv"
    if not os.path.isdir(foldercoord):
        os.makedirs(foldercoord)
        print("Folder created:", name)
        print("saved csv:", name)
    else:

        print("saved csv:", name)
    df.to_csv(os.path.join (foldercoord, name))






from genericpath import exists

from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import umap
import mdtraj as mdt
import numpy as np
import glob as glob
import pandas as pd
import os as os
import umap
import matplotlib
import matplotlib.pyplot as plt
from math import sqrt
from matplotlib.lines import Line2D
from pyemma.coordinates import tica

from matplotlib.lines import Line2D
import cv2
import imageio.v2 as imageio
import shutil




def getnums(s, e,i):
   return list(range(s, e,i))


def dihedrals_angles(root,proteinpath):

    for pt in proteinpath:
        print(os.path.join(root,pt,"*.xtc"))
        #idangles = {"idCA":[], "idHN":[], "idN":[],"idHA":[],"idresi":[],"nameresi":[]}
        tjs=glob.glob(os.path.join(root,pt,"*.xtc"))
        print(tjs)
        for xtc in tjs:
            print(xtc)
            pdb = glob.glob(os.path.join(root,pt,"*.pdb"))[0] # Geometric estructure of the
            print(pdb)
            # Read in xtc file
            traj = mdt.load(xtc, top=pdb)
            traj=traj.remove_solvent()
            with open(pdb) as f:
                lines = f.readlines()

            lst= getnums(50,2501,50)
            fr=0
            for nframemax in lst:
                idangles={"idCA":[], "idHN":[], "idN":[],"idHA":[],"idresi":[],"nameresi":[],"idresidues":[]}
                #print(nframemax)
                subtraj=traj[int(nframemax-50):int(nframemax)]
                for frame in subtraj:
                    t = frame.topology
                    for i in t.atoms:
                        if i.name == "CA":
                            #print(i.residue)
                            idangles["idCA"].append(i.index)
                            idangles["idHN"].append(i.index-1)
                            idangles["idN"].append(i.index-2)
                            idangles["idHA"].append(i.index+1)
                            idangles["idresi"].append(str(i.residue)[3:])
                            idangles["nameresi"].append(str(i.residue)[:3])
                            idangles["idresidues"].append(int(fr))
                            fr=fr+1

                df=pd.DataFrame(idangles)
                print(df.shape)
                indices=np.array(df)
                indices=indices[:,:4] #
                diedrals= mdt.compute_dihedrals(subtraj, indices, periodic=True, opt=True) # in radianes
                #pd.DataFrame(diedrals)
                psi=mdt.compute_psi(subtraj, periodic=True, opt=True)
                phi=mdt.compute_phi(subtraj, periodic=True, opt=True)
                chi=mdt.compute_chi1(subtraj, periodic=True, opt=True)
                folderresrmsd = os.path.join(root,pt, "coordinates", "angles")
                namedied=xtc.split("/")[-1].split(".xtc")[0]+"_"+ str(nframemax)+"_died.csv"
                namepsi1=xtc.split("/")[-1].split(".xtc")[0]+"_"+ str(nframemax)+"_index_psi.csv"
                namepsi2=xtc.split("/")[-1].split(".xtc")[0]+"_"+ str(nframemax)+"_apsi.csv"
                namephi1=xtc.split("/")[-1].split(".xtc")[0]+"_"+ str(nframemax)+"_index_phi.csv"
                namephi2=xtc.split("/")[-1].split(".xtc")[0]+"_"+ str(nframemax)+"_aphi.csv"
                namechi1=xtc.split("/")[-1].split(".xtc")[0]+"_"+ str(nframemax)+"_index_chi.csv"
                namechi2=xtc.split("/")[-1].split(".xtc")[0]+"_"+ str(nframemax)+"_achi.csv"

                if not os.path.isdir(folderresrmsd):
                    os.makedirs(folderresrmsd)
                    print("created")
                else:
                    print("exist:")

                pd.DataFrame(diedrals).to_csv(os.path.join(root,pt,"coordinates", "angles",namedied))
                pd.DataFrame(psi[0]).to_csv(os.path.join(root,pt,"coordinates", "angles",namepsi1))
                pd.DataFrame(psi[1]).to_csv(os.path.join(root,pt,"coordinates", "angles",namepsi2))
                pd.DataFrame(phi[0]).to_csv(os.path.join(root,pt,"coordinates", "angles",namephi1))
                pd.DataFrame(phi[1]).to_csv(os.path.join(root,pt,"coordinates", "angles",namephi2))
                pd.DataFrame(chi[0]).to_csv(os.path.join(root,pt,"coordinates", "angles",namechi1))
                pd.DataFrame(chi[1]).to_csv(os.path.join(root,pt,"coordinates", "angles",namechi2))
                #TODO crear carpeta para los csvs. y elimininar los que hay

def resspoints(root,proteinpath):
    for pt in proteinpath:
        trajlst=glob.glob(os.path.join(root,pt,"coordinates","*res.csv"))

        dic={"trayectory":[],"protein":[],"framesid":[],"idresidue":[],"massx":[],"massy":[],"massz":[]}
        for traj in trajlst[:]:
            print(traj)
            df=pd.read_csv(traj)
            #print(df.head)
            for row in df.itertuples():
                dic["trayectory"].append(traj.split("/")[-1].split("_")[0])
                dic["protein"].append(traj.split("_")[-3].split("res.csv")[0])
                #dic["frames"].append(traj.split("_")[-1].split("res.csv")[0])
                dic["framesid"].append(row.idframe)
                #print(row.idframe)
                dic["idresidue"].append(row.idress)
                dic["massx"].append(row.massx)
                dic["massy"].append(row.massy)
                dic["massz"].append(row.massz)

        dfcoord=pd.DataFrame(dic)
        dfcoord = dfcoord.dropna()
        folderresrmsd = os.path.join(root,pt, "coordinates", "3Drespositions")
        name="3Dres_"+pt+".csv"
        if not os.path.isdir(folderresrmsd):
            os.makedirs(folderresrmsd)
            print("saved csv:", name)
        else:
            print("saved csv:", name)
        dfcoord.to_csv(os.path.join(folderresrmsd, name))
        #dfcoord.sort_values(["idresidue"], ascending = True)

def helix (root, proteinlst, lsthelixs):
  for protein in proteinlst[:]:
        pt=protein
        namecsv="3Dres_"+pt+".csv"
        dfh=pd.read_csv(os.path.join(root,pt, "coordinates", "3Drespositions",namecsv))
        dfh = dfh.sort_values(by = 'framesid')
        print(dfh.iloc[:,3:])
        #dfh = dfcoord[(dfcoord['idresidue'].astype(int) == int(res))]#| (dfcoord['idresidue'].astype(int) == int(32)) | (dfcoord['idresidue'].astype(int) == int(60)) ] #| dfcoord['idresidue'].astype(int) == int(ah1-1) |dfcoord['idresidue'].astype(int) == int(ah1+1) ]
        for trajec in dfh.trayectory.unique().tolist()[:]:
            print(trajec)
            dicrmsd={"framesid":[],"idposition":[],"idresidue1":[],"idresidue2":[],"protein":[],"RMSD":[],"trayectory":[]}
            dfhtraj = dfh[(dfh["trayectory"].astype(int) == int(trajec))]
            #print(dfhtraj)
            for frame in dfhtraj.framesid.unique().tolist()[:]:
                #print(frame)
                dfhf = dfhtraj[(dfhtraj['framesid'].astype(int) == int(frame))]
                #print(dfhf)
                for res in range(3)[:]:
                    for helix in range(7)[:]:
                        print(helix,res)
                        print(lsthelixs[helix][res])
                        dfhi1=dfhf[(dfhf['idresidue'].astype(int) == int(lsthelixs[helix][res]))]
                        print("SDFafsd")
                        for i in range(7):
                            helix=0
                            dfhi2=dfhf[(dfhf['idresidue'].astype(int) == int(lsthelixs[helix+i][res]))]
                            #print(lsthelixs[helix+i])

                            #print(dfhi1.idresidue.values[0],dfhi2.idresidue.values[0])

                            RMSD=sqrt((dfhi2.massx.values[0]-dfhi1.massx.values[0]) ** 2)+ \
                                     ((dfhi2.massy.values[0]-dfhi1.massy.values[0]) ** 2) + \
                                     ((dfhi2.massz.values[0]-dfhi1.massz.values[0]) ** 2)
                            dicrmsd["framesid"].append(frame)
                            #print(frame)

                            dicrmsd["idposition"].append(res)
                            dicrmsd["idresidue1"].append(int(dfhi1.idresidue.values[0]))
                            dicrmsd["idresidue2"].append(int(dfhi2.idresidue.values[0]))
                            dicrmsd["protein"].append(pt)
                            dicrmsd["trayectory"].append(dfhi2.trayectory.values[0])
                            dicrmsd["RMSD"].append(RMSD)
                            #print(dfhi1.idresidue.values[0],dfhi2.idresidue.values[0])

            df=pd.DataFrame(dicrmsd)
            #df4= df[((df.idresidue1 == int(46)) &( df.idresidue2 == int(46) ))].index
            #df=df.drop(df[df.idresidue1 == df.idresidue2].index)

            dr= df[((df.idresidue1 == int(32)) &( df.idresidue2 == int(70)))| \
                   ((df.idresidue1 == int(32)) &( df.idresidue2 == int(106) )) | \
                    ((df.idresidue1 == int(32)) &( df.idresidue2 == int(150) )) | \
                    ((df.idresidue1 == int(32)) &( df.idresidue2 == int(200) )) | \
                    ((df.idresidue1 == int(32)) &( df.idresidue2 == int(284) )) | \
                   ((df.idresidue1 == int(32)) &( df.idresidue2 == int(322) ))  | \

                    ((df.idresidue1 == int(70)) &( df.idresidue2 == int(106) )) | \
                    ((df.idresidue1 == int(70)) &( df.idresidue2 == int(150) )) | \
                    ((df.idresidue1 == int(70)) &( df.idresidue2 == int(200) )) | \
                    ((df.idresidue1 == int(70)) &( df.idresidue2 == int(284) )) | \
                    ((df.idresidue1 == int(70)) &( df.idresidue2 == int(322) )) | \

                    ((df.idresidue1 == int(106)) &( df.idresidue2 == int(150) )) | \
                     ((df.idresidue1 == int(106)) &( df.idresidue2 == int(200) )) | \
                    ((df.idresidue1 == int(106)) &( df.idresidue2 == int(284) )) | \
                    ((df.idresidue1 == int(106)) &( df.idresidue2 == int(322) )) | \

                     ((df.idresidue1 == int(150)) &( df.idresidue2 == int(200) )) | \
                     ((df.idresidue1 == int(150)) &( df.idresidue2 == int(284) )) | \
                    ((df.idresidue1 == int(150)) &( df.idresidue2 == int(322) )) | \

                     ((df.idresidue1 == int(200)) &( df.idresidue2 == int(284) )) | \
                     ((df.idresidue1 == int(200)) &( df.idresidue2 == int(322) )) | \

                    ((df.idresidue1 == int(284)) &( df.idresidue2 == int(322) )) | \

                #####################################################################

                    ((df.idresidue1 == int(46)) &( df.idresidue2 == int(83)))| \
                   ((df.idresidue1 == int(46)) &( df.idresidue2 == int(121) )) | \
                    ((df.idresidue1 == int(46)) &( df.idresidue2 == int(160) )) | \
                    ((df.idresidue1 == int(46)) &( df.idresidue2 == int(215) )) | \
                    ((df.idresidue1 == int(46)) &( df.idresidue2 == int(298) )) | \
                   ((df.idresidue1 == int(46)) &( df.idresidue2 == int(331) ))  | \

                    ((df.idresidue1 == int(83)) &( df.idresidue2 == int(121) )) | \
                    ((df.idresidue1 == int(83)) &( df.idresidue2 == int(160) )) | \
                    ((df.idresidue1 == int(83)) &( df.idresidue2 == int(215) )) | \
                    ((df.idresidue1 == int(83)) &( df.idresidue2 == int(298) )) | \
                    ((df.idresidue1 == int(83)) &( df.idresidue2 == int(331) )) | \

                    ((df.idresidue1 == int(121)) &( df.idresidue2 == int(160) )) | \
                     ((df.idresidue1 == int(121)) &( df.idresidue2 == int(215) )) | \
                    ((df.idresidue1 == int(121)) &( df.idresidue2 == int(298) )) | \
                    ((df.idresidue1 == int(121)) &( df.idresidue2 == int(331) )) | \

                     ((df.idresidue1 == int(160)) &( df.idresidue2 == int(215) )) | \
                     ((df.idresidue1 == int(160)) &( df.idresidue2 == int(298) )) | \
                    ((df.idresidue1 == int(160)) &( df.idresidue2 == int(331) )) | \

                     ((df.idresidue1 == int(215)) &( df.idresidue2 == int(298) )) | \
                     ((df.idresidue1 == int(215)) &( df.idresidue2 == int(331) )) | \

                    ((df.idresidue1 == int(298)) &( df.idresidue2 == int(331) )) | \

                #####################################################################

                    ((df.idresidue1 == int(60)) &( df.idresidue2 == int(96)))| \
                   ((df.idresidue1 == int(60)) &( df.idresidue2 == int(136) )) | \
                    ((df.idresidue1 == int(60)) &( df.idresidue2 == int(171) )) | \
                    ((df.idresidue1 == int(60)) &( df.idresidue2 == int(229) )) | \
                    ((df.idresidue1 == int(60)) &( df.idresidue2 == int(312) )) | \
                   ((df.idresidue1 == int(60)) &( df.idresidue2 == int(341) ))  | \

                    ((df.idresidue1 == int(96)) &( df.idresidue2 == int(136) )) | \
                    ((df.idresidue1 == int(96)) &( df.idresidue2 == int(171) )) | \
                    ((df.idresidue1 == int(96)) &( df.idresidue2 == int(229) )) | \
                    ((df.idresidue1 == int(96)) &( df.idresidue2 == int(312) )) | \
                    ((df.idresidue1 == int(96)) &( df.idresidue2 == int(341) )) | \

                    ((df.idresidue1 == int(136)) &( df.idresidue2 == int(171) )) | \
                     ((df.idresidue1 == int(136)) &( df.idresidue2 == int(229) )) | \
                    ((df.idresidue1 == int(136)) &( df.idresidue2 == int(312) )) | \
                    ((df.idresidue1 == int(136)) &( df.idresidue2 == int(341) )) | \

                     ((df.idresidue1 == int(171)) &( df.idresidue2 == int(229) )) | \
                     ((df.idresidue1 == int(171)) &( df.idresidue2 == int(312) )) | \
                    ((df.idresidue1 == int(171)) &( df.idresidue2 == int(341) )) | \

                     ((df.idresidue1 == int(229)) &( df.idresidue2 == int(312) )) | \
                     ((df.idresidue1 == int(229)) &( df.idresidue2 == int(341) )) | \

                    ((df.idresidue1 == int(312)) &( df.idresidue2 == int(341) ))  ]


            name=protein+"_"+str(trajec) +"_helix.csv"
            namecenter=protein+"_"+str(trajec) +"_centerhelix.csv"

            folder=os.path.join(root,pt,"coordinates","helix")
            dc= dr[dr.idposition == int(1)]

            if not os.path.isdir(folder):
                os.makedirs(folder)
                print("saved csv:", folder)
            else:
                print("exist:", folder)
            dr.to_csv(os.path.join(folder, name))
            dc.to_csv(os.path.join(folder, namecenter))

def extract_angles(datapath,prot):
    dic={"trayectory":[], "protein":[], "idframe":[], "idresidue":[], "idangle":[], "valueangle":[]}
    for file in glob.glob(datapath)[:]:
        print(file)
        if not "index" in file:
            if not "died" in file:
                print(file)
                #df = pd.read_csv(file)
                #print(df.shape)
                if not "achi" in file:
                    df = pd.read_csv(file)
                    df=df.iloc[:,1:]#TODO camboiar sacar 30
                    #print(df.shape)
                    name=file.split("/")[-1]
                    trajname=name[:5]
                    for res in range(df.shape[1]):
                        fr=0
                        for j in range (df.shape[0])[:]:#TODO camboiar sacar 10
                            dic["trayectory"].append(trajname)
                            dic["protein"].append(name.split("_")[2])
                            dic["idangle"].append(name.split("_")[-1].split(".")[0])
                            framenum=int(name.split("_")[3]) -50 + fr
                            dic["idframe"].append(framenum)
                            print(framenum)
                            fr=fr+1
                            dic["idresidue"].append(res)
                            ang=df.iloc[j][res]
                            dic["valueangle"].append((ang))
    df=pd.DataFrame(dic)
    #y = "{:.2f}".format(x)
    folder=os.path.join(datapath.split("*.csv")[0],"embeddings")
    if not os.path.isdir(folder):
        os.makedirs(folder)
        print("saved csv:", folder)
    else:
        print("saved csv:", folder)
    df.to_csv(os.path.join(folder,"angles_"+prot+".csv"))

def UMAP_angles(datapath,prot):
    root=datapath.split("*.csv")[0]
    name="angles_"+prot+".csv"
    dfprot=pd.read_csv(os.path.join(root,"embeddings",name))
    dfprot=dfprot.sort_values(by=["idresidue",'idframe'])

    test_list=dfprot.idframe.unique().tolist()
    #len(test_list)
    #for i, j in zip(test_list,range(0,2500)):
       # print(i,j)

    print("vericar que esten todos los frames",len(test_list))
    for traj in dfprot.trayectory.unique().tolist()[:]:
        print(traj)
        df1= dfprot[dfprot.trayectory.astype(int) == int(traj)]
        #dfaphi.iloc[:,3:]
        dfapsi= df1[df1.idangle.astype(str) == "apsi"]
        dfaphi= df1[df1.idangle.astype(str) == "aphi"]

        rowspsi=[]
        for fr in dfapsi.idframe.unique().tolist()[:]:
            dfframepsi= dfapsi[dfapsi.idframe.astype(int) == int(fr)]
            dfframepsi = dfframepsi.set_index('idresidue').T
            datapsi=rowspsi.append(dfframepsi.iloc[5:,0:])
        dfpsi = pd.concat(rowspsi)
        #dicnew={"apsi":dfapsi.valueangle.tolist(),"aphi": dfaphi.valueangle.tolist()}
        #dfn=pd.DataFrame(dicnew)
        # Creating DataFrame by passing Dictionary
        #data = dfn[["apsi", "aphi"]].values
        data=dfpsi.values
        reducer = umap.UMAP()
        scaled_data = StandardScaler().fit_transform(data)
        print(scaled_data.shape)
        embedding = reducer.fit_transform(scaled_data)
        print(embedding.shape)
        nam= str(traj)+'_psi_umap_allresidues.csv'
        folder=os.path.join(root,"embeddings","allresidues","UMAP")
        if not os.path.isdir(folder):
            os.makedirs(folder)
            print("saved csv:", folder)
        else:
            print("saved csv:", folder)
        np.savetxt(os.path.join(folder,nam), embedding, delimiter=",")
    ##################################################
        rowsphi=[]
        for fr in dfaphi.idframe.unique().tolist()[:]:
            dfframephi= dfaphi[dfaphi.idframe.astype(int) == int(fr)]
            dfframephi = dfframephi.set_index('idresidue').T
            dataphi=rowsphi.append(dfframephi.iloc[5:,0:])
        dfphi = pd.concat(rowsphi)
        #dicnew={"apsi":dfapsi.valueangle.tolist(),"aphi": dfaphi.valueangle.tolist()}
        #dfn=pd.DataFrame(dicnew)
        # Creating DataFrame by passing Dictionary
        #data = dfn[["apsi", "aphi"]].values
        data=dfphi.values
        reducer = umap.UMAP()
        scaled_data = StandardScaler().fit_transform(data)
        print(scaled_data.shape)
        embedding = reducer.fit_transform(scaled_data)
        print(embedding.shape)
        nam= str(traj)+'_phi_umap_allresidues.csv'
        np.savetxt(os.path.join(folder,nam), embedding, delimiter=",")

        ######################
        dfbothangle=pd.concat([dfphi,dfpsi],axis=1)
        data=dfbothangle.values
        reducer = umap.UMAP()
        scaled_data = StandardScaler().fit_transform(data)
        print(scaled_data.shape)
        embedding = reducer.fit_transform(scaled_data)
        print(embedding.shape)
        nam= str(traj)+'_psi_phi_umap_allresidues.csv'
        np.savetxt(os.path.join(folder,nam), embedding, delimiter=",")


def UMAP_helix(datapath):

    for traj in glob.glob(datapath)[:]:
        print(traj)
        if not "center" in traj:
            if not "Icon"in traj:
                print(traj)
                df=pd.read_csv(traj)
                print(df)
                dfpos=[]
                for position in df.idposition.unique().tolist()[:]:
                    print(position)
                    df1= df[df.idposition.astype(int) == int(position)]
                    df1=df1.sort_values(by=["idresidue1",'framesid'])
                    #df1.iloc[:,:4]
                    rows=[]
                    for fr in df1.framesid.unique().tolist()[:]:#TODO verificar porque hay 2400
                        df2= df1[df1.framesid.astype(int) == int(fr)]
                        df2 = df2.set_index('idresidue1').T
                        rows.append(df2.iloc[5:6,0:])
                    df3 = pd.concat(rows)
                    dfpos.append(df3.values)
                    data=df3.values
                    reducer = umap.UMAP()
                    scaled_data = StandardScaler().fit_transform(data)
                    print(scaled_data.shape)
                    embedding = reducer.fit_transform(scaled_data)
                    print(embedding.shape)
                    folder=os.path.join(datapath.split("*.csv")[0],"embeddings","UMAP")
                    if not os.path.isdir(folder):
                        os.makedirs(folder)
                        print("saved csv:", folder)
                    else:
                        print("saved csv:", folder)
                    nam=os.path.join(folder,str(traj.split("/")[-1].split("_")[-2])+"_"+str(position)+'_helix_umap.csv')
                    print(nam)
                    np.savetxt(nam, embedding, delimiter=",")
                ##############all the positions###########
                data = np.concatenate((dfpos[0],dfpos[1],dfpos[2]),axis=1)
                reducer = umap.UMAP()
                scaled_data = StandardScaler().fit_transform(data)
                print(scaled_data.shape)
                embedding = reducer.fit_transform(scaled_data)
                print(embedding.shape)
                nam=os.path.join(folder,str(traj.split("/")[-1].split("_")[-2])+"_"+str("allpositions")+'_helix_umap.csv')
                print(nam)
                np.savetxt(nam, embedding, delimiter=",")
def UMAP_3Dposition(lstrprot):
    for prot in lstrprot[:]:
        datapath="./dataset/"+prot+"/coordinates/3Drespositions/3Dres_"+str(prot)+".csv"
        df=pd.read_csv(datapath)
        df["idresidue"]=df["idresidue"].astype(int)
        df=df.sort_values(by=["idresidue",'framesid'])


        for traj in df.trayectory.unique().tolist()[:]:

            print(traj)
            df1= df[df.trayectory.astype(int) == int(traj)]
            print(df1)

            rowsx=[]
            rowsy=[]
            rowsz=[]
            for fr in df1.framesid.unique().tolist()[:]:
                df2= df1[df1.framesid.astype(int) == int(fr)]
                df2 = df2.set_index('idresidue').T
                rowsx.append(df2.iloc[4:5,0:])
                rowsy.append(df2.iloc[5:6,0:])
                rowsz.append(df2.iloc[6:7,0:])
            dfx = pd.concat(rowsx)
            dfy = pd.concat(rowsy)
            dfz = pd.concat(rowsz)



            datax=dfx.values
            datay=dfy.values
            dataz=dfy.values

            reducerx = umap.UMAP()
            reducery = umap.UMAP()
            reducerz = umap.UMAP()

            scaled_datax = StandardScaler().fit_transform(datax)
            scaled_datay = StandardScaler().fit_transform(datay)
            scaled_dataz = StandardScaler().fit_transform(dataz)

            print(scaled_datax.shape)
            print(scaled_datay.shape)
            print(scaled_dataz.shape)

            embeddingx = reducerx.fit_transform(scaled_datax)
            embeddingy = reducery.fit_transform(scaled_datay)
            embeddingz = reducerz.fit_transform(scaled_dataz)

            print(embeddingx.shape)
            print(embeddingy.shape)
            print(embeddingz.shape)

            folder=os.path.join(datapath.split("3Dres_"+str(prot)+".csv")[0],"embeddings", "allresidues","UMAP")
            if not os.path.isdir(folder):
                os.makedirs(folder)
                print("saved csv:", folder)
            else:
                print("saved csv:", folder)

            namx=os.path.join(folder,str(traj)+'_3Dx_umap.csv')
            namy=os.path.join(folder,str(traj)+'_3Dy_umap.csv')
            namz=os.path.join(folder,str(traj)+'_3Dz_umap.csv')

            np.savetxt(namx, embeddingx, delimiter=",")
            np.savetxt(namy, embeddingy, delimiter=",")
            np.savetxt(namz, embeddingz, delimiter=",")
            ##############with xyz positions###########

            dataxzy = np.concatenate((dfx.values,dfy.values,dfz.values),axis=1)

            scaled_dataxyz = StandardScaler().fit_transform(dataxzy)
            print(scaled_dataxyz.shape)
            reducerxzy = umap.UMAP()
            embeddingxyz = reducerxzy.fit_transform(scaled_dataxyz)
            print(embeddingxyz.shape)
            namxyz=os.path.join(folder,str(traj)+'_3Dxyz_umap.csv')
            np.savetxt(namxyz, embeddingxyz, delimiter=",")
            print(namxyz)


def TICA_3Dposition(lstrprot):
    for prot in lstrprot[:]:
            datapath="./dataset/"+prot+"/coordinates/3Drespositions/3Dres_"+str(prot)+".csv"
            df=pd.read_csv(datapath)
            df["idresidue"]=df["idresidue"].astype(int)
            df=df.sort_values(by=["idresidue",'framesid'])


            for traj in df.trayectory.unique().tolist()[:]:

                print(traj)
                df1= df[df.trayectory.astype(int) == int(traj)]
                #print(df1)

                rowsx=[]
                rowsy=[]
                rowsz=[]
                for fr in df1.framesid.unique().tolist()[:]:
                    df2= df1[df1.framesid.astype(int) == int(fr)]
                    df2 = df2.set_index('idresidue').T
                    rowsx.append(df2.iloc[4:5,0:])
                    rowsy.append(df2.iloc[5:6,0:])
                    rowsz.append(df2.iloc[6:7,0:])
                dfx = pd.concat(rowsx)
                dfy = pd.concat(rowsy)
                dfz = pd.concat(rowsz)



                datax=dfx.values
                datay=dfy.values
                dataz=dfy.values


                print(datax.shape)
                print(datay.shape)
                print(dataz.shape)

                embeddingx = tica(datax, lag=2, dim=2).get_output()[0]
                embeddingy = tica(datay, lag=2, dim=2).get_output()[0]
                embeddingz = tica(dataz, lag=2, dim=2).get_output()[0]



                print(embeddingx.shape)
                print(embeddingy.shape)
                print(embeddingz.shape)

                folder=os.path.join(datapath.split("3Dres_"+str(prot)+".csv")[0],"embeddings", "allresidues","TICA")
                if not os.path.isdir(folder):
                    os.makedirs(folder)
                    print("saved csv:", folder)
                else:
                    print("saved csv:", folder)

                namx=os.path.join(folder,str(traj)+'_3Dx_tica.csv')
                namy=os.path.join(folder,str(traj)+'_3Dy_tica.csv')
                namz=os.path.join(folder,str(traj)+'_3Dz_tica.csv')

                np.savetxt(namx, embeddingx, delimiter=",")
                np.savetxt(namy, embeddingy, delimiter=",")
                np.savetxt(namz, embeddingz, delimiter=",")
                ##############with xyz positions###########

                dataxzy = np.concatenate((dfx.values,dfy.values,dfz.values),axis=1)


                embeddingxyz = tica(dataxzy, lag=2, dim=2).get_output()[0]
                print(embeddingxyz.shape)
                namxyz=os.path.join(folder,str(traj)+'_3Dxyz_tica.csv')
                np.savetxt(namxyz, embeddingxyz, delimiter=",")
                print(namxyz)
def TICA_angles(lstrprot):
    for prot in lstrprot[:]:
        datapath="./dataset/"+prot+"/coordinates/angles/*.csv"
        root=datapath.split("*.csv")[0]
        name="angles_"+prot+".csv"
        dfprot=pd.read_csv(os.path.join(root,"embeddings",name))
        dfprot=dfprot.sort_values(by=["idresidue",'idframe'])

        test_list=dfprot.idframe.unique().tolist()
        #len(test_list)
        #for i, j in zip(test_list,range(0,2500)):
           # print(i,j)

        print("vericar que esten todos los frames",len(test_list))
        for traj in dfprot.trayectory.unique().tolist()[:]:
            print(traj)
            df1= dfprot[dfprot.trayectory.astype(int) == int(traj)]
            #dfaphi.iloc[:,3:]
            dfapsi= df1[df1.idangle.astype(str) == "apsi"]
            dfaphi= df1[df1.idangle.astype(str) == "aphi"]

            rowspsi=[]
            for fr in dfapsi.idframe.unique().tolist()[:]:
                dfframepsi= dfapsi[dfapsi.idframe.astype(int) == int(fr)]
                dfframepsi = dfframepsi.set_index('idresidue').T
                datapsi=rowspsi.append(dfframepsi.iloc[5:,0:])
            dfpsi = pd.concat(rowspsi)
            #dicnew={"apsi":dfapsi.valueangle.tolist(),"aphi": dfaphi.valueangle.tolist()}
            #dfn=pd.DataFrame(dicnew)
            # Creating DataFrame by passing Dictionary
            #data = dfn[["apsi", "aphi"]].values
            data=dfpsi.values
            data=np.float64(data)
            embedding = tica(data, lag=2, dim=2).get_output()[0]
            nam= str(traj)+'_psi_tica_allresidues.csv'
            folder=os.path.join(root,"embeddings","allresidues","TICA")
            if not os.path.isdir(folder):
                os.makedirs(folder)
                print("saved csv:", folder)
            else:
                print("saved csv:", folder)
            np.savetxt(os.path.join(folder,nam), embedding, delimiter=",")
        ##################################################
            rowsphi=[]
            for fr in dfaphi.idframe.unique().tolist()[:]:
                dfframephi= dfaphi[dfaphi.idframe.astype(int) == int(fr)]
                dfframephi = dfframephi.set_index('idresidue').T
                dataphi=rowsphi.append(dfframephi.iloc[5:,0:])
            dfphi = pd.concat(rowsphi)
            #dicnew={"apsi":dfapsi.valueangle.tolist(),"aphi": dfaphi.valueangle.tolist()}
            #dfn=pd.DataFrame(dicnew)
            # Creating DataFrame by passing Dictionary
            #data = dfn[["apsi", "aphi"]].values
            data=dfphi.values
            data=np.float64(data)
            embedding = tica(data, lag=2, dim=2).get_output()[0]
            nam= str(traj)+'_phi_tica_allresidues.csv'
            np.savetxt(os.path.join(folder,nam), embedding, delimiter=",")
            ######################
            dfbothangle=pd.concat([dfphi,dfpsi],axis=1)
            data=dfbothangle.values
            data=np.float64(data)
            embedding = tica(data, lag=2, dim=2).get_output()[0]
            nam= str(traj)+'_psi_phi_tica_allresidues.csv'
            np.savetxt(os.path.join(folder,nam), embedding, delimiter=",")

def TICA_helix(datapath):

    for traj in glob.glob(datapath)[:]:
        print(traj)
        if not "center" in traj:
            if not "Icon"in traj:
                print(traj)
                df=pd.read_csv(traj)
                print(df)
                dfpos=[]
                for position in df.idposition.unique().tolist()[:]:
                    print(position)
                    df1= df[df.idposition.astype(int) == int(position)]
                    df1=df1.sort_values(by=["idresidue1",'framesid'])
                    #df1.iloc[:,:4]
                    rows=[]
                    for fr in df1.framesid.unique().tolist()[:]:#TODO verificar porque hay 2400
                        df2= df1[df1.framesid.astype(int) == int(fr)]
                        df2 = df2.set_index('idresidue1').T
                        rows.append(df2.iloc[5:6,0:])
                    df3 = pd.concat(rows)
                    dfpos.append(df3.values)
                    data=df3.values
                    data=np.array(data, dtype=float)
                    print(data.shape)
                    data=np.float64(data)
                    embedding = tica(data, lag=2, dim=2).get_output()[0]
                    print(embedding.shape)
                    folder=os.path.join(datapath.split("*.csv")[0],"embeddings","TICA")
                    if not os.path.isdir(folder):
                        os.makedirs(folder)
                        print("saved csv:", folder)
                    else:
                        print("saved csv:", folder)
                    nam=os.path.join(folder,str(traj.split("/")[-1].split("_")[-2])+"_"+str(position)+'_helix_tica.csv')
                    print(nam)
                    np.savetxt(nam, embedding, delimiter=",")
                ##############all the positions###########
                data = np.concatenate((dfpos[0],dfpos[1],dfpos[2]),axis=1)
                data=np.array(data, dtype=float)
                print(data.shape)
                data=np.float64(data)
                embedding = tica(data, lag=2, dim=2).get_output()[0]
                print(embedding.shape)
                nam=os.path.join(folder,str(traj.split("/")[-1].split("_")[-2])+"_"+str("allpositions")+'_helix_tica.csv')
                print(nam)
                np.savetxt(nam, embedding, delimiter=",")
def graph (datapath):
    for path in glob.glob(datapath)[:]:
        name = path.split("/")[-1].split(".csv")[0]
        if not "angles" in name:
            #print(path)
            #df = pd.read_csv(path)
            df = pd.read_csv(path, header=None)
            pathindex = os.path.join((path.split("embeddings")[0]),"caroline_distances","RMSD", path.split("/")[-1].split("_")[0]+"RMSD.csv")
            #pathindex = pathindex.replace("angles", "3Drespositions")
            print(pathindex)
            dfindex = pd.read_csv(pathindex)
            print(dfindex)

            #sizes = [20*2**n  for n in range(df.shape[0])]
            #sizes = [(i*0.01)+1 for i in range(df.shape[0])]

            #df = df[dfindex.state.astype(str) == "inactive"]

            #dfa = df[dfindex.state.astype(str) == "active"]
            #dft = df[dfindex.state.astype(str) == "transition"]

            x1 = df.iloc[:, 0]
            y1 = df.iloc[:, 1]
            size = [(i*0.05)+1 for i in range(x1.shape[0])]

            x2 = dfindex.iloc[:, 0]
            y2 = dfindex.iloc[:, 2]

            #plt.figure()
            #plt.figure().clear()
            #plt.figure()

            col = []
            st = []
            dfv = dfindex.state.values
            raw = dfindex.RMSD.values.tolist()
            amin, amax = min(raw), max(raw)

            for i, val in enumerate(raw):
                raw[i] = ((val-amin) / (amax-amin))

            for i in range(0, len(dfv)):
                if  dfv[i] == "active":
                    col.append('#7E2F8E')
                    st.append('active')
                elif dfv[i] == "inactive":
                    col.append('#EDB120')
                    st.append('inactive')
                elif dfv[i] == "transition":
                    col.append('#D95319')
                    st.append('transition')

            # Crear figura y subplots
            fig, axs = plt.subplots(2)
            # Ajustar tamaño de la figura
            fig.set_size_inches(10, 7)
            fig.suptitle(name, fontweight='heavy', fontsize=20, x=0.5)
            fig.tight_layout(pad=2, h_pad=4)

            # Primer gráfico
            axs[0].scatter(x1, y1, color=col, s=size, alpha=raw, edgecolors='black', linewidth=0.3)
            legend_elements = [
                Line2D([0], [0], marker='o', color='w', label='Active', markerfacecolor='#7E2F8E', markersize=8),
                Line2D([0], [0], marker='o', color='w', label='Inactive', markerfacecolor='#EDB120', markersize=8),
                Line2D([0], [0], marker='o', color='w', label='Transition', markerfacecolor='#D95319', markersize=8)
                ]
            axs[0].legend(handles=legend_elements, loc='best')
            axs[0].set_title('Título del subplot 1', fontweight='heavy', loc='center')
            axs[0].set_xlabel('Eje X')
            axs[0].set_ylabel('Eje Y')
            axs[0].spines['top'].set_visible(False)
            axs[0].spines['right'].set_visible(False)
            axs[0].grid(color='grey', alpha=0.2)
            # Segundo gráfico
            axs[1].plot(x2, y2, linewidth=1, label='RMSD')
            axs[1].grid(color='grey', alpha=0.2)
            axs[1].legend(loc='best')
            axs[1].set_title('Título del subplot 2', fontweight='heavy', loc='center')
            axs[1].set_xlabel('Columna 1')
            axs[1].set_ylabel('Columna 2')

            #scatter plot
            #for i in range(len(x)):
                #plt.scatter(x[i],y[i], color=col[i],s=size[i], edgecolors='black', linewidth=0.3)
            #plt.gca().set_aspect('equal', 'datalim')

            plt.savefig(path.split(".csv")[0]+".png")
            print(path.split(".csv")[0]+".png")
            plt.close()
def graph_line(datapath):
    for path in glob.glob(datapath)[:]:
        print(path)
        #print(glob.glob(datapath)[:3])
        name = path.split("/")[-1].split(".csv")[0]
        if not "angles" in name:
            print(path)
            #df = pd.read_csv(path)
            df = pd.read_csv(path, header=None)
            pathindex = os.path.join((path.split("embeddings")[0]),"caroline_distances","RMSD", path.split("/")[-1].split("_")[0]+"RMSD.csv")
            #pathindex = pathindex.replace("angles", "3Drespositions")
            print(pathindex)
            if pathindex.split("/")[-1] != "RMSD.csv":

                dfindex = pd.read_csv(pathindex)
                print(dfindex)

                #sizes = [20*2**n  for n in range(df.shape[0])]
                #sizes = [(i*0.01)+1 for i in range(df.shape[0])]

                #df = df[dfindex.state.astype(str) == "inactive"]

                #dfa = df[dfindex.state.astype(str) == "active"]pyplot_show
                #dft = df[dfindex.state.astype(str) == "transition"]

                x1 = df.iloc[:, 0]
                y1 = df.iloc[:, 1]
                size = [(i*0.05)+1 for i in range(x1.shape[0])]

                x2 = dfindex.iloc[:, 0]
                y2 = dfindex.iloc[:, 2]

                #plt.figure()
                #plt.figure().clear()
                #plt.figure()

                col = []
                st = []
                dfv = dfindex.state.values
                raw = dfindex.RMSD.values.tolist()
                amin, amax = min(raw), max(raw)

                for i, val in enumerate(raw):
                    raw[i] = ((val-amin) / (amax-amin))

                for i in range(0, len(dfv)):
                    if  dfv[i] == "active":
                        col.append('#7E2F8E')
                        st.append('active')
                    elif dfv[i] == "inactive":
                        col.append('#EDB120')
                        st.append('inactive')
                    elif dfv[i] == "transition":
                        col.append('#D95319')
                        st.append('transition')

                # Crear figura y subplots
                fig, axs = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [4, 2]})
                fig.suptitle(str(name), fontweight='heavy', fontsize=20, x=0.5) #name cambiar el name
                fig.tight_layout(pad=2, h_pad=4)

                # Primer gráfico
                axs[0].set_aspect('equal')
                axs[0].scatter(x1, y1, color=col, s=size, alpha=raw, edgecolors='black', linewidth=0.3)
                axs[0].grid(color='grey', alpha=0.2)
                legend_elements = [
                    Line2D([0], [0], marker='o', color='w', label='Active', markerfacecolor='#7E2F8E', markersize=8),
                    Line2D([0], [0], marker='o', color='w', label='Inactive', markerfacecolor='#EDB120', markersize=8),
                    Line2D([0], [0], marker='o', color='w', label='Transition', markerfacecolor='#D95319', markersize=8)
                    ]
                axs[0].legend(handles=legend_elements, loc='best')
                a=[axs[0].set_xlim(),axs[0].set_ylim()]
                values = [val for tup in a for val in tup]
                axs[0].set_xlim((min(values),max(values)))
                axs[0].set_title('Título del subplot 1', fontweight='heavy', loc='center')
                axs[0].set_xlabel('Eje X')
                axs[0].set_ylabel('Eje Y')
                axs[0].spines['top'].set_visible(False)
                axs[0].spines['right'].set_visible(False)

                # Segundo gráficox
                mask_active = dfindex['state'] == 'active'
                mask_inactive = dfindex['state'] == 'inactive'
                mask_transition = dfindex['state'] == 'transition'
                col2 = {'active': '#7E2F8E', 'inactive': '#EDB120', 'transition': '#D95319'}
                for i in range(len(x2) - 1):
                    xi, yi = x2[i], y2[i]
                    xf, yf = x2[i+1], y2[i+1]
                    if mask_active[i]:
                        seg = 'active'
                    elif mask_inactive[i]:
                        seg = 'inactive'
                    elif mask_transition[i]:
                        seg = 'transition'
                    axs[1].plot([xi, xf], [yi, yf], color=col2[seg])
                axs[1].grid(color='grey', alpha=0.2)
                legend_elements2 = [
                    Line2D([0], [0], color='#7E2F8E', label='Active', markersize=8),
                    Line2D([0], [0], color='#EDB120', label='Inactive', markersize=8),
                    Line2D([0], [0], color='#D95319', label='Transition', markersize=8)
                    ]
                axs[1].legend(handles=legend_elements2, loc='best')
                axs[1].set_title('Título del subplot 2', fontweight='heavy', loc='center')
                axs[1].set_xlabel('Columna 1')
                axs[1].set_ylabel('Columna 2')

                #scatter plot
                #for i in range(len(x)):
                    #plt.scatter(x[i],y[i], color=col[i],s=size[i], edgecolors='black', linewidth=0.3)
                #plt.gca().set_aspect('equal', 'datalim')

                plt.savefig(path.split(".csv")[0]+".png")
                print(path.split(".csv")[0]+".png")
                plt.close()
            else:
                continue
def graph_area(datapath):
    for path in glob.glob(datapath)[:]:
        name = path.split("/")[-1].split(".csv")[0]
        if not "angles" in name:
            #print(path)
            #df = pd.read_csv(path)
            df = pd.read_csv(path, header=None)
            pathindex = os.path.join((path.split("embeddings")[0]),"caroline_distances","RMSD", path.split("/")[-1].split("_")[0]+"RMSD.csv")
            #pathindex = pathindex.replace("angles", "3Drespositions")
            print(pathindex)


            dfindex = pd.read_csv(pathindex)
            print(dfindex)

            #sizes = [20*2**n  for n in range(df.shape[0])]
            #sizes = [(i*0.01)+1 for i in range(df.shape[0])]

            #df = df[dfindex.state.astype(str) == "inactive"]

            #dfa = df[dfindex.state.astype(str) == "active"]
            #dft = df[dfindex.state.astype(str) == "transition"]

            x1 = df.iloc[:, 0]
            y1 = df.iloc[:, 1]
            size = [(i*0.05)+1 for i in range(x1.shape[0])]

            x2 = dfindex.iloc[:, 0]
            y2 = dfindex.iloc[:, 2]

            #plt.figure()
            #plt.figure().clear()
            #plt.figure()

            col = []
            st = []
            dfv = dfindex.state.values
            raw = dfindex.RMSD.values.tolist()
            amin, amax = min(raw), max(raw)

            for i, val in enumerate(raw):
                raw[i] = ((val-amin) / (amax-amin))

            for i in range(0, len(dfv)):
                if  dfv[i] == "active":
                    col.append('#7E2F8E')
                    st.append('active')
                elif dfv[i] == "inactive":
                    col.append('#EDB120')
                    st.append('inactive')
                elif dfv[i] == "transition":
                    col.append('#D95319')
                    st.append('transition')

            # Crear figura y subplots
            fig, axs = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [4, 2]})
            fig.suptitle(name, fontweight='heavy', fontsize=20, x=0.5)
            fig.tight_layout(pad=2, h_pad=4)

             # Primer gráfico
            axs[0].set_aspect('equal')
            axs[0].scatter(x1, y1, color=col, s=size, alpha=raw, edgecolors='black', linewidth=0.3)
            axs[0].grid(color='grey', alpha=0.2)
            legend_elements = [
                Line2D([0], [0], marker='o', color='w', label='Active', markerfacecolor='#7E2F8E', markersize=8),
                Line2D([0], [0], marker='o', color='w', label='Inactive', markerfacecolor='#EDB120', markersize=8),
                Line2D([0], [0], marker='o', color='w', label='Transition', markerfacecolor='#D95319', markersize=8)
                ]
            axs[0].legend(handles=legend_elements, loc='best')
            a=[axs[0].set_xlim(),axs[0].set_ylim()]
            values = [val for tup in a for val in tup]
            axs[0].set_xlim((min(values),max(values)))
            axs[0].set_title('Título del subplot 1', fontweight='heavy', loc='center')
            axs[0].set_xlabel('Eje X')
            axs[0].set_ylabel('Eje Y')
            axs[0].spines['top'].set_visible(False)
            axs[0].spines['right'].set_visible(False)

            # Segundo gráfico
            state_colors = {'active': '#7E2F8E', 'inactive': '#EDB120', 'transition': '#D95319'}
            for (xi, yi), (xf, yf), state in zip(zip(x2, y2), zip(x2[1:], y2[1:]), dfindex['state']):
                axs[1].plot([xi, xf], [yi, yf], color=state_colors[state])
                # Colorear a partir de un gráfico de área
                axs[1].fill_between([xi, xf], [yi, yf], color=state_colors[state], alpha=0.2)
            axs[1].grid(color='grey', alpha=0.2)
            legend_elements2 = [
                Line2D([0], [0], color='#7E2F8E', label='Active', markersize=8),
                Line2D([0], [0], color='#EDB120', label='Inactive', markersize=8),
                Line2D([0], [0], color='#D95319', label='Transition', markersize=8)
                ]
            axs[1].legend(handles=legend_elements2, loc='best')
            axs[1].set_title('Título del subplot 2', fontweight='heavy', loc='center')
            axs[1].set_xlabel('Columna 1')
            axs[1].set_ylabel('Columna 2')

            #scatter plot
            #for i in range(len(x)):
                #plt.scatter(x[i],y[i], color=col[i],s=size[i], edgecolors='black', linewidth=0.3)
            #plt.gca().set_aspect('equal', 'datalim')
            name=path.split(".csv")[0]+".png"

            rute=path.split("/")[:-1]
            folder = os.path.join(*rute,"area")
            if not os.path.isdir(folder):
                os.makedirs(folder)
                print("saved csv:", folder)
            else:
                print("saved csv:", folder)

            rute=path.split("/")[:-1]
            nam=name.split("/")[-1]

            plt.savefig(os.path.join(folder,nam))

            plt.close()
        else:
            print(pathindex)
def graph_gift(datapath):
    for path in glob.glob(datapath)[:]:
        print(path)
        name = path.split("/")[-1].split(".csv")[0]
        if not "angles" in name:
            df = pd.read_csv(path, header=None)
            print(path)
            pathindex = os.path.join((path.split("embeddings")[0]), "caroline_distances", "RMSD", path.split("/")[-1].split("_")[0] + "RMSD.csv")
            print(pathindex)

            if pathindex.split("/")[-1] != "RMSD.csv":
                dfindex = pd.read_csv(pathindex)
                print(dfindex)

                x1, y1 = df.iloc[:, 0].to_numpy(), df.iloc[:, 1].to_numpy()
                size = [(i * 0.05) + 1 for i in range(x1.shape[0])]
                x2, y2 = dfindex.iloc[:, 0].to_numpy(), dfindex.iloc[:, 2].to_numpy()


                dfv = dfindex.state.values
                col = np.empty(len(dfv), dtype=object)
                raw = dfindex['RMSD'].to_numpy()
                amin, amax = raw.min(), raw.max()
                raw = (raw - amin) / (amax - amin)

                col[dfv == "active"] = '#7E2F8E'
                col[dfv == "inactive"] = '#EDB120'
                col[dfv == "transition"] = '#D95319'

                legend_elements = [
                    Line2D([0], [0], marker='o', color='w', label='Active', markerfacecolor='#7E2F8E', markersize=8),
                    Line2D([0], [0], marker='o', color='w', label='Inactive', markerfacecolor='#EDB120', markersize=8),
                    Line2D([0], [0], marker='o', color='w', label='Transition', markerfacecolor='#D95319', markersize=8)
                ]

                legend_elements2 = [
                    Line2D([0], [0], color='#7E2F8E', label='Active', markersize=8),
                    Line2D([0], [0], color='#EDB120', label='Inactive', markersize=8),
                    Line2D([0], [0], color='#D95319', label='Transition', markersize=8)
                ]

                # Crear figura, subplots y configuración
                fig, axs = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [4, 2]})
                fig.suptitle(str(name), fontweight='heavy', fontsize=20, x=0.5)
                fig.tight_layout(pad=2, h_pad=4)

                # Primer gráfico
                axs[0].scatter(x1, y1, color=col, s=size, alpha=raw, edgecolors='black', linewidth=0.3)
                axs[0].set_aspect('equal')
                #longest_axis = max(max(x1) - min(x1), max(y2) - min(y2))
                #if longest_axis is x1:
                    #axs[0].set_ylim(axs[0].get_xlim())
                #else:
                    #axs[0].set_xlim(axs[0].get_ylim())

                axs[0].grid(color='grey', alpha=0.2)
                #a=[axs[0].set_xlim(),axs[0].set_ylim()]
                #values = [val for tup in a for val in tup]
                axs[0].set_xlim(min(x1),max(x1))
                axs[0].set_ylim(min(y1),max(y1))
                axs[0].legend(handles=legend_elements, loc='best')
                axs[0].set_title('Posiciones X e Y para cada amonoácido', fontweight='heavy', loc='center')
                axs[0].set_xlabel('Eje X')
                axs[0].set_ylabel('Eje Y')
                axs[0].spines['top'].set_visible(False)
                axs[0].spines['right'].set_visible(False)

                # Segundo gráfico
                state_colors = {'active': '#7E2F8E', 'inactive': '#EDB120', 'transition': '#D95319'}
                for (xi, yi), (xf, yf), state in zip(zip(x2, y2), zip(x2[1:], y2[1:]), dfindex['state']):
                    axs[1].plot([xi, xf], [yi, yf], color=state_colors[state])
                    axs[1].fill_between([xi, xf], [yi, yf], color=state_colors[state], alpha=0.04)
                axs[1].grid(color='grey', alpha=0.2)
                axs[1].legend(handles=legend_elements2, loc='best')
                axs[1].set_title('Distancias', fontweight='heavy', loc='center')
                axs[1].set_xlabel('Tiempo (frames)')
                axs[1].set_ylabel('RMSD')

                # Guardar figura estática
                plt.savefig(path.split(".csv")[0]+".png", transparent=False, facecolor='white')
                plt.close()

                # Crear carpeta para guardar las imágenes para el video
                folder = os.path.join(os.path.dirname(path), str(name))
                if not os.path.isdir(folder):
                    os.makedirs(folder)
                    print("Carpeta creada exitosamente.")
                else:
                    print("La carpeta ya existe.")

                def generar_imgs():
                # Graficar dinámicamente los datos
                    def graf_din(n, num_points):
                        fig, axs = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [4, 2]})
                        fig.suptitle(str(name), fontweight='heavy', fontsize=20, x=0.5)
                        fig.tight_layout(pad=2, h_pad=4)

                        axs[0].scatter(x1[:num_points], y1[:num_points], color=col[:num_points], s=size[:num_points], alpha=raw[:num_points], edgecolors='black', linewidth=0.3)
                        axs[0].set_aspect('equal')
                        axs[0].set_xlim(min(x1)-int(5),max(x1)+int(5))
                        axs[0].set_ylim(min(y1)-int(5),max(y1)+int(5))
                        axs[0].grid(color='grey', alpha=0.2)
                        axs[0].legend(handles=legend_elements, loc='best')
                        axs[0].set_title('Dimensionality reduction', fontweight='heavy', loc='center')
                        axs[0].set_xlabel('Eje X')
                        axs[0].set_ylabel('Eje Y')
                        axs[0].spines['top'].set_visible(False)
                        axs[0].spines['right'].set_visible(False)

                        for (xi, yi), (xf, yf), state in zip(
                            zip(x2[:num_points], y2[:num_points]),
                            zip(x2[1:num_points + 1], y2[1:num_points + 1]),
                            dfindex['state'][:num_points]):
                                axs[1].plot([xi, xf], [yi, yf], color=state_colors[state])
                                axs[1].fill_between([xi, xf], [yi, yf], color=state_colors[state], alpha=0.04)
                        axs[1].grid(color='grey', alpha=0.2)
                        axs[1].legend(handles=legend_elements2, loc='best')
                        axs[1].set_title('Distancias', fontweight='heavy', loc='center')
                        axs[1].set_xlabel('Tiempo (frames)')
                        axs[1].set_ylabel('RMSD')

                        plt.savefig(os.path.join(folder, name + f'_img_{n}.png'), transparent=False, facecolor='white')
                        print(os.path.join(folder, name + f'_img_{n}.png'))
                        plt.close()

                    num_points = 50
                    for i in range(num_points, len(x1) + 1, num_points):
                        graf_din(i, i)

                generar_imgs()

                # Crear video y gif
                images = sorted([img for img in os.listdir(folder) if img.endswith(".png")], key=lambda x: int(x.split("_")[-1].split(".")[0]))

                frame = cv2.imread(os.path.join(folder, images[0]))
                height, width, _ = frame.shape

                fourcc = cv2.VideoWriter_fourcc(*'avc1')
                videofolder=path.split("/")[:-1]
                videofolder=os.path.join(*videofolder,"videos")
                if not os.path.isdir(videofolder):
                    os.makedirs(videofolder)
                    print("Carpeta creada exitosamente.")
                else:
                    print("La carpeta ya existe.")
                namevideo=path.split("/")[-1].split(".csv")[0]
                pathvideo=os.path.join(videofolder,namevideo)
                video_lento = cv2.VideoWriter(pathvideo + "_lento.mp4", fourcc, 2, (width, height))
                #video_media = cv2.VideoWriter(pathvideo + "_media.mp4", fourcc, 6, (width, height))
                #video_rapido = cv2.VideoWriter(pathvideo + "_rapido.mp4", fourcc, 15, (width, height))

                frames = []
                for image in images:
                    image_path = os.path.join(folder, image)
                    frame = cv2.imread(image_path)
                    video_lento.write(frame)
                    #video_media.write(frame)
                    #video_rapido.write(frame)
                    frames.append(imageio.imread(image_path))

                cv2.destroyAllWindows()
                video_lento.release()
                #video_media.release()
                #video_rapido.release()

                imageio.mimsave(pathvideo + ".gif", frames, 'GIF', loop=0)

                # Eliminar las imágenes luego de crear el video y el gif
                shutil.rmtree(folder)
            else:
                continue
def PCA(X , num_components):

    #Step-1
    X_meaned = X - np.mean(X , axis = 0)

    #Step-2
    cov_mat = np.cov(X_meaned , rowvar = False)

    #Step-3
    eigen_values , eigen_vectors = np.linalg.eigh(cov_mat)

    #Step-4
    sorted_index = np.argsort(eigen_values)[::-1]
    sorted_eigenvalue = eigen_values[sorted_index]
    sorted_eigenvectors = eigen_vectors[:,sorted_index]

    #Step-5
    eigenvector_subset = sorted_eigenvectors[:,0:num_components]

    #Step-6
    X_reduced = np.dot(eigenvector_subset.transpose() , X_meaned.transpose() ).transpose()

    return X_reduced
def PCA_angles(datapath,prot):
    root=datapath.split("*.csv")[0]
    name="angles_"+prot+".csv"
    dfprot=pd.read_csv(os.path.join(root,"embeddings",name))
    dfprot=dfprot.sort_values(by=["idresidue",'idframe'])

    test_list=dfprot.idframe.unique().tolist()
    #len(test_list)
    #for i, j in zip(test_list,range(0,2500)):
       # print(i,j)

    print("vericar que esten todos los frames",len(test_list))
    for traj in dfprot.trayectory.unique().tolist()[:]:
        print(traj)
        df1= dfprot[dfprot.trayectory.astype(int) == int(traj)]
        #dfaphi.iloc[:,3:]
        dfapsi= df1[df1.idangle.astype(str) == "apsi"]
        dfaphi= df1[df1.idangle.astype(str) == "aphi"]

        rowspsi=[]
        for fr in dfapsi.idframe.unique().tolist()[:]:
            dfframepsi= dfapsi[dfapsi.idframe.astype(int) == int(fr)]
            dfframepsi = dfframepsi.set_index('idresidue').T
            rowspsi.append(dfframepsi.iloc[5:,0:])
        dfpsi = pd.concat(rowspsi)
        #dicnew={"apsi":dfapsi.valueangle.tolist(),"aphi": dfaphi.valueangle.tolist()}
        #dfn=pd.DataFrame(dicnew)
        # Creating DataFrame by passing Dictionary
        #data = dfn[["apsi", "aphi"]].values
        data=dfpsi.values
        data=np.array(data, dtype=float)
        print(data.shape)
        a=PCA(data,num_components=2) # YOU CAN  CHOSE THE COMPONENT NUMBER
        print("psi: ",a.shape)
        #nam= str(traj)+'_psi_PCA_allresidues.csv'
        folder=os.path.join(root,"embeddings","allresidues","PCA")
        if not os.path.isdir(folder):
            os.makedirs(folder)
            print("saved csv:", folder)
        else:
            print("saved csv:", folder)
        nam=os.path.join(folder,str(traj)+'_psi_PCA_allresidues.csv')
        np.savetxt(nam, a, delimiter=",")
    ##################################################
        rowsphi=[]
        for fr in dfaphi.idframe.unique().tolist()[:]:
            dfframephi= dfaphi[dfaphi.idframe.astype(int) == int(fr)]
            dfframephi = dfframephi.set_index('idresidue').T
            rowsphi.append(dfframephi.iloc[5:,0:])
        dfphi = pd.concat(rowsphi)
        #dicnew={"apsi":dfapsi.valueangle.tolist(),"aphi": dfaphi.valueangle.tolist()}
        #dfn=pd.DataFrame(dicnew)
        # Creating DataFrame by passing Dictionary
        #data = dfn[["apsi", "aphi"]].values
        data=dfphi.values
        data=np.array(data, dtype=float)
        print(data.shape)
        p=PCA(data,num_components=2) # YOU CAN  CHOSE THE COMPONENT NUMBER
        print("p: ",p.shape)
        nam=os.path.join(folder,str(traj)+'_phi_PCA_allresidues.csv')
        np.savetxt(nam, p, delimiter=",")
        ######################
        dfbothangle=pd.concat([dfphi,dfpsi],axis=1)
        data=dfbothangle.values
        data=np.array(data, dtype=float)
        print(data.shape)
        phipsi=PCA(data,num_components=2) # YOU CAN  CHOSE THE COMPONENT NUMBER
        print(phipsi.shape)
        nam=os.path.join(folder,str(traj)+'_phi_psi_PCA_allresidues.csv')
        np.savetxt(nam, phipsi, delimiter=",")
def PCA_helix(datapath):

    for traj in glob.glob(datapath)[:]:
        print(traj)
        if not "center" in traj:
            if not "Icon"in traj:
                print(traj)
                df=pd.read_csv(traj)
                print(df)
                dfpos=[]
                for position in df.idposition.unique().tolist()[:]:
                    print(position)
                    df1= df[df.idposition.astype(int) == int(position)]
                    df1=df1.sort_values(by=["idresidue1",'framesid'])
                    #df1.iloc[:,:4]
                    rows=[]
                    for fr in df1.framesid.unique().tolist()[:]:#TODO verificar porque hay 2400
                        df2= df1[df1.framesid.astype(int) == int(fr)]
                        df2 = df2.set_index('idresidue1').T
                        rows.append(df2.iloc[5:6,0:])
                    df3 = pd.concat(rows)
                    dfpos.append(df3.values)
                    data=df3.values
                    data=np.array(data, dtype=float)
                    print(data.shape)
                    a=PCA(data,num_components=2) # YOU CAN  CHOSE THE COMPONENT NUMBER
                    print("psi: ",a.shape)
                    folder=os.path.join(datapath.split("*.csv")[0],"embeddings","PCA")
                    if not os.path.isdir(folder):
                        os.makedirs(folder)
                        print("saved csv:", folder)
                    else:
                        print("saved csv:", folder)
                    nam=os.path.join(folder,str(traj.split("/")[-1].split("_")[-2])+"_"+str(position)+'_helix_PCA.csv')
                    print(nam)
                    np.savetxt(nam, a, delimiter=",")
                ##############all the positions###########
                data = np.concatenate((dfpos[0],dfpos[1],dfpos[2]),axis=1)
                data=np.array(data, dtype=float)
                print(data.shape)
                b=PCA(data,num_components=2) # YOU CAN  CHOSE THE COMPONENT NUMBER
                print("helix: ",b.shape)
                nam=os.path.join(folder,str(traj.split("/")[-1].split("_")[-2])+"_"+str("allpositions")+'_helix_PCA.csv')
                print(nam)
                np.savetxt(nam, b, delimiter=",")
def PCA_3Dposition(lstrprot):
    for prot in lstrprot[:]:
        datapath="./dataset/"+prot+"/coordinates/3Drespositions/3Dres_"+str(prot)+".csv"
        df=pd.read_csv(datapath)
        df["idresidue"]=df["idresidue"].astype(int)
        df=df.sort_values(by=["idresidue",'framesid'])


        for traj in df.trayectory.unique().tolist()[:]:
            print(traj)
            df1= df[df.trayectory.astype(int) == int(traj)]
            print(df1)

            rowsx=[]
            rowsy=[]
            rowsz=[]
            for fr in df1.framesid.unique().tolist()[:]:
                df2= df1[df1.framesid.astype(int) == int(fr)]
                df2 = df2.set_index('idresidue').T
                rowsx.append(df2.iloc[4:5,0:])
                rowsy.append(df2.iloc[5:6,0:])
                rowsz.append(df2.iloc[6:7,0:])
            dfx = pd.concat(rowsx)
            dfy = pd.concat(rowsy)
            dfz = pd.concat(rowsz)



            datax=dfx.values
            datay=dfy.values
            dataz=dfy.values

            datax=np.array(datax, dtype=float)
            datay=np.array(datay, dtype=float)
            dataz=np.array(dataz, dtype=float)

            print(datax.shape)
            print(datay.shape)
            print(dataz.shape)


            ax=PCA(datax,num_components=2) # YOU CAN  CHOSE THE COMPONENT NUMBER
            ay=PCA(datay,num_components=2) # YOU CAN  CHOSE THE COMPONENT NUMBER
            az=PCA(dataz,num_components=2) # YOU CAN  CHOSE THE COMPONENT NUMBER

            print("ax: ",ax.shape)
            print("ay: ",ay.shape)
            print("az: ",az.shape)



            folder=os.path.join(datapath.split("3Dres_"+str(prot)+".csv")[0],"embeddings","allresidues","PCA")
            if not os.path.isdir(folder):
                os.makedirs(folder)
                print("saved csv:", folder)
            else:
                print("saved csv:", folder)

            namx=os.path.join(folder,str(traj)+'_3Dx_PCA.csv')
            namy=os.path.join(folder,str(traj)+'_3Dy_PCA.csv')
            namz=os.path.join(folder,str(traj)+'_3Dz_PCA.csv')

            np.savetxt(namx, ax, delimiter=",")
            np.savetxt(namy, ay, delimiter=",")
            np.savetxt(namz, az, delimiter=",")
            ##############with xyz positions###########

            dataxzy = np.concatenate((dfx.values,dfy.values,dfz.values),axis=1)

            dataxzy=np.array(dataxzy, dtype=float)
            print("dataxzy: ",dataxzy.shape)
            axyz=PCA(dataxzy,num_components=2) # YOU CAN  CHOSE THE COMPONENT NUMBER
            namxyz=os.path.join(folder,str(traj)+'_3Dxyz_PCA.csv')
            print("dataxzy: ",axyz.shape)
            np.savetxt(namxyz, axyz, delimiter=",")
def PCA_3Dpositionresidues(lstrprot):
    for prot in lstrprot[:]:
        datapath="./dataset/"+prot+"/coordinates/3Drespositions/3Dres_"+str(prot)+".csv"
        df=pd.read_csv(datapath)
        df["idresidue"]=df["idresidue"].astype(int)
        df=df.sort_values(by=["idresidue",'framesid'])
        for traj in df.trayectory.unique().tolist()[:]:
            print(traj)
            df1= df[(df.trayectory.astype(int) == int(traj))& (df.idresidue.astype(int) <= int(272)) & (df.idresidue.astype(int) >= int(131))]
            print(df1)

            rowsx=[]
            rowsy=[]
            rowsz=[]
            for fr in df1.framesid.unique().tolist()[:]:
                df2= df1[df1.framesid.astype(int) == int(fr)]
                df2 = df2.set_index('idresidue').T
                rowsx.append(df2.iloc[4:5,0:])
                rowsy.append(df2.iloc[5:6,0:])
                rowsz.append(df2.iloc[6:7,0:])
            dfx = pd.concat(rowsx)
            dfy = pd.concat(rowsy)
            dfz = pd.concat(rowsz)



            datax=dfx.values
            datay=dfy.values
            dataz=dfy.values

            datax=np.array(datax, dtype=float)
            datay=np.array(datay, dtype=float)
            dataz=np.array(dataz, dtype=float)

            print(datax.shape)
            print(datay.shape)
            print(dataz.shape)


            ax=PCA(datax,num_components=2) # YOU CAN  CHOSE THE COMPONENT NUMBER
            ay=PCA(datay,num_components=2) # YOU CAN  CHOSE THE COMPONENT NUMBER
            az=PCA(dataz,num_components=2) # YOU CAN  CHOSE THE COMPONENT NUMBER

            print("ax: ",ax.shape)
            print("ay: ",ay.shape)
            print("az: ",az.shape)



            folder=os.path.join(datapath.split("3Dres_"+str(prot)+".csv")[0],"embeddings","residues", "PCA")
            if not os.path.isdir(folder):
                os.makedirs(folder)
                print("saved csv:", folder)
            else:
                print("saved csv:", folder)

            namx=os.path.join(folder,str(traj)+'_3Dresiduesx_PCA.csv')
            namy=os.path.join(folder,str(traj)+'_3Dresiduesy_PCA.csv')
            namz=os.path.join(folder,str(traj)+'_3Dresiduesz_PCA.csv')

            np.savetxt(namx, ax, delimiter=",")
            np.savetxt(namy, ay, delimiter=",")
            np.savetxt(namz, az, delimiter=",")
            ##############with xyz positions###########

            dataxzy = np.concatenate((dfx.values,dfy.values,dfz.values),axis=1)

            dataxzy=np.array(dataxzy, dtype=float)
            print("dataxzy: ",dataxzy.shape)
            axyz=PCA(dataxzy,num_components=2) # YOU CAN  CHOSE THE COMPONENT NUMBER
            namxyz=os.path.join(folder,str(traj)+'_3Dresiduesxyz_PCA.csv')
            print("dataxzy: ",axyz.shape)
            np.savetxt(namxyz, axyz, delimiter=",")


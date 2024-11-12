

datapath="./dataset/*/coordinates/3Drespositions/embeddings/*/*/*.csv"
fn.graph(datapath)

def graph (datapath):
    for path in glob.glob(datapath)[:]:
        name=path.split("/")[-1].split(".csv")[0]
        if not "angles" in name:
            #print(path)
            #df=pd.read_csv(path)
            df = pd.read_csv(path, header=None)
            pathindex=os.path.join((path.split("embeddings")[0]),"distances",path.split("/")[-1].split("_")[0]+".csv")
            #pathindex=pathindex.replace("angles", "3Drespositions")
            print(pathindex)
            dfindex = pd.read_csv(pathindex)
            print(dfindex)


            #sizes = [20*2**n  for n in range(df.shape[0])]
            #sizes = [(i*0.01)+1 for i in range(df.shape[0])]

            #df= df[dfindex.state.astype(str) == "inactive"]

            #dfa= df[dfindex.state.astype(str) == "active"]
            #dft= df[dfindex.state.astype(str) == "transition"]

            x = df.iloc[:,0]
            y = df.iloc[:,1]
            #plt.figure()
            #plt.figure().clear()
            #plt.figure()
            size = [(i*0.05)+1 for i in range(x.shape[0])]

            #color list with same length as the data
            col=[]
            st=[]
            dfv=dfindex.state.values
            raw=dfindex.RMSD.values.tolist()
            amin, amax = min(raw), max(raw)
            for i, val in enumerate(raw):
                raw[i] = ((val-amin) / (amax-amin))

            for i in range(0,len(dfv)):
                if  dfv[i] == "active":
                    col.append('#7E2F8E')
                    st.append('active')
                elif dfv[i] == "inactive":
                    col.append('#EDB120')
                    st.append('inactive')
                elif dfv[i] == "transition":
                    col.append('#D95319')
                    st.append('transition')

            plt.scatter(x,y, color=col ,s=size, alpha=raw, edgecolors='black', linewidth=0.3)

            legend_elements = [Line2D([0], [0], marker='o', color='w', label='Active',
                                      markerfacecolor='#7E2F8E', markersize=8),
                               Line2D([0], [0], marker='o', color='w', label='Inactive',
                                      markerfacecolor='#EDB120', markersize=8),
                               Line2D([0], [0], marker='o', color='w', label='Transition',
                                      markerfacecolor='#D95319', markersize=8)]

            plt.legend(handles=legend_elements, loc='best')

            #scatter plot
            #for i in range(len(x)):
                #plt.scatter(x[i],y[i], color=col[i],s=size[i], edgecolors='black', linewidth=0.3)
            #plt.gca().set_aspect('equal', 'datalim')
            plt.title(name)
            plt.xlabel('Columna 1')
            plt.ylabel('Columna 2')
            plt.savefig(path.split(".csv")[0]+".png")
            print(path.split(".csv")[0]+".png")
            plt.close()

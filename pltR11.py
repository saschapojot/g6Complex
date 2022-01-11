import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


#load shooting data
prefix="startL4"
shootingDf=pd.read_csv(prefix+"shootingR11.csv")

gShooting=shootingDf["g"]
EShooting=shootingDf["E"]

###loading upper data
levelNum=1
upperDf=pd.read_csv("level"+str(levelNum)+"upperGComplexR11.csv")
gUpper=upperDf["g"]
EReUpper=upperDf["ERe"]
EImUpper=upperDf["EIm"]



################
EMax=80
indLarge=np.logical_and(np.abs(EImUpper)<EMax,np.abs(EImUpper)>1)
indSmall=np.where(np.abs(EImUpper)<1)[0]
gUpperLarge=gUpper[indLarge]
gUpperSmall=gUpper[indSmall]
EReUpperLarge=EReUpper[indLarge]
EReUpperSmall=EReUpper[indSmall]
EImUpperLarge=EImUpper[indLarge]
EImUpperSmall=EImUpper[indSmall]

########
fig, ax = plt.subplots(figsize=(20, 20))
ax.set_ylabel("E")
# plt.yscale('symlog')
ax.set_xscale("log")
ax.set_xlabel("g")
ax.set_title("Eigenvalues for potential $V(x)=ix^{3}-gx^{4}$")


sctShooting=ax.scatter(gShooting,EShooting,color="blue",marker=".",s=40,label="Shooting")
sctUpperReSmall=ax.scatter(gUpperSmall,EReUpperSmall,color="fuchsia",marker="+",s=50,label="WKB upper real |im|<1")
sctUpperImSmall=ax.scatter(gUpperSmall,EImUpperSmall,color="green",marker="^",s=50,label="WKB upper imag |im|<1")
sctUpperReLarge=ax.scatter(gUpperLarge,EReUpperLarge,color="brown",marker="+",s=50,label="WKB upper real |im|>1")
sctUpperImLarge=ax.scatter(gUpperLarge,EImUpperLarge,color="cyan",marker="^",s=50,label="WKB upper imag |im|>1")
ax.hlines(y=0,xmin=1e-4,xmax=10,color="r")
plt.legend()
plt.savefig("n="+str(levelNum)+"tmp11.png")
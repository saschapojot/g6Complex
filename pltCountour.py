import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#this script plots wkb results using contour integral


#load shooting data
prefix="num"+str(100)+"startL10"
shootingDf=pd.read_csv(prefix+"deleted0.1shootingR11.csv")

gShooting=shootingDf["g"]
EShooting=shootingDf["E"]

#wkb data
inWKBFile="./level0end7contourUpperGComplexR11.csv"


upperDf=pd.read_csv(inWKBFile)

gUpper=np.array(upperDf["g"])
EReUpper=np.array(upperDf["ERe"])
EImUpper=np.array(upperDf["EIm"])

######
EMax=80
indReSmall=np.where(np.abs(EReUpper)<EMax)

gSmall=gUpper[indReSmall]
EReSmall=EReUpper[indReSmall]
EImSmall=EImUpper[indReSmall]

########
fig, ax = plt.subplots(figsize=(20, 20))
ax.set_ylabel("E")
# plt.yscale('symlog')
ax.set_xscale("log")
ax.set_xlabel("g")
ax.set_title("Eigenvalues for potential $V(x)=ix^{3}-gx^{4}$")



sctShooting=ax.scatter(gShooting,EShooting,color="blue",marker=".",s=40,label="Shooting")
sctSmall=ax.scatter(gSmall,EReSmall,color="fuchsia",marker="+",label="WKB",s=40)
plt.legend()
plt.savefig("contour.png")
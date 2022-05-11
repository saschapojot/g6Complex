from funcsContour import *
import pandas as pd
from multiprocessing import Pool
from datetime import datetime


num=5
startG=1e-4

stopG=1e0
gnIndAll = mpmath.linspace(mpmath.log10(startG), mpmath.log10(stopG), num)

gAll = [10 ** elem for elem in gnIndAll]
threadNum = 48
levelStart=0
levelEnd=2
levelsAll = range(levelStart, levelEnd + 1)
inDataAll=[]

for nTmp in levelsAll:
    for gTmp in gAll:
        EEst=(nTmp+1/2)*mp.pi
        inDataAll.append([nTmp,gTmp,EEst])


############################################################################
# ###########parallel computation  for upper, may be memory consuming
tWKBUpperStart=datetime.now()
pool1=Pool(threadNum)
retAllUpper=pool1.map(computeOneSolutionWithUpperPairs,inDataAll)
tWKBAdjEnd=datetime.now()
print("parallel WKB time for upper pairs: ",tWKBAdjEnd-tWKBUpperStart)

######data serialization for upper
nValsUpper=[]
gValsUpper=[]
ERealValsUpper=[]
EImagValsUpper=[]
for itemTmp in retAllUpper:
    if len(itemTmp)==0:
        continue
    n,g,E=itemTmp
    # delete complex solutions
    if np.abs(np.imag(E)) > 1e-3:
        continue
    EReTmp=mpmath.re(E)
    EImTmp=mpmath.im(E)
    nValsUpper.append(n)
    gValsUpper.append(g)
    ERealValsUpper.append(EReTmp)
    EImagValsUpper.append(EImTmp)

###write data of upper to csv
upperDat=np.array([nValsUpper,gValsUpper,ERealValsUpper,EImagValsUpper]).T
upperDf=pd.DataFrame(upperDat,columns=["n","g","ERe","EIm"])
upperDf.to_csv("level"+str(levelStart)+"end"+str(levelEnd)+"contourUpperGComplexR11.csv",index=False)
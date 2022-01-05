from shootingR11Funcs import *
import pandas as pd
import scipy.special as sspecial

####values of g and E
num = 500
startG = 1e-4
stopG = 1e1

gnIndAll = np.linspace(start=np.log10(startG), stop=np.log10(stopG), num=num)
gAll = [10 ** elem for elem in gnIndAll]
EMax = 40.0
# convert to lambda and F
inDataAll = []  # contains [lambda, FEst]

dE = 0.1
for g in gAll:
    for E in np.arange(0.1, EMax, dE):
        lmd = g ** (-5 / 6)
        F = E * g ** (-1 / 3)
        inDataAll.append([lmd, F])

threadNum = 24
pool1 = Pool(threadNum)
tShootingStart = datetime.now()
retAll = pool1.map(computeOneSolution, inDataAll)
tShootingEnd = datetime.now()
print("Shooting time: ", tShootingEnd - tShootingStart)

# data serialization
gShootingVals = []
EShootingVals = []
for itemTmp in retAll:
    lmd, F = itemTmp
    if F < 0:
        continue
    gTmp = lmd ** (-6 / 5)
    ETmp = F * lmd ** (-2 / 5)
    if np.abs(ETmp) > 100:
        continue
    gShootingVals.append(gTmp)
    EShootingVals.append(ETmp)

##plot shooting
fig, ax = plt.subplots(figsize=(20, 20))
ax.set_xscale("log")
ax.set_ylabel("E")
ax.set_xlabel("g")
ax.set_title("Shooting eigenvalues for potential $V(x)=ix^{3}-gx^{4}$, region I-II")

shootingScatter = ax.scatter(gShootingVals, EShootingVals, color="blue", marker=".", s=50, label="shooting")

# lowerWKB5Sct=ax.scatter(gWKBSct,EWKBSct,color="red",marker="x",s=40,label="WKB lower $-igx^{5}$")
plt.legend()

plt.savefig("startL" + str(LEst) + "shootingR12" + str(dE) + "deleted.png")
plt.close()

dataPdFrame = np.array([gShootingVals, EShootingVals]).T
dfgE = pd.DataFrame(dataPdFrame, columns=["g", "E"])
dfgE.to_csv("startL" + str(LEst) + "shootingR12.csv")

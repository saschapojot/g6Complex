from datetime import  datetime
from multiprocessing import Pool
import mpmath
from mpmath import  mp
import  numpy as np
mp.dps=12

#this script computes potential 6 V=ix^{3}-gx^{4}

def retUpperPair(g,E):
    """

    :param g:
    :param E:
    :return: the pair of roots whose im > 0, the first root has re <0, the second
    has re >0
    """
    coef=[-np.float(g),1j,0,0,-np.complex(E)]

    # rootsAll=mpmath.polyroots(coef,maxsteps=100)
    rootsAll=np.roots(coef)
    tmp=[]
    for elem in rootsAll:
        if mpmath.im(elem)>0:
            tmp.append(elem)

    tmp=sorted(tmp,key=mpmath.re)
    rst=[tmp]
    return rst



def f(z,g,E):
    """

    :param z:
    :param g:
    :param E:
    :return:(-Q)^{1/2}
    """
    return mpmath.sqrt(g*z**4-1j*z**3+E)



def integralQuadrature(g,E,x1,x2):
    """

    :param g:
    :param E:
    :param x1: ending point
    :param x2: starting point
    :return:
    """
    a1 = mpmath.re(x1)
    b1 = mpmath.im(x1)

    a2 = mpmath.re(x2)
    b2 = mpmath.im(x2)
    slope = (b1 - b2) / (a1 - a2)
    gFunc = lambda y: f(y + 1j * (slope * (y - a2) + b2), g, E)
    return (1 + 1j * slope) * mpmath.quad(gFunc, [a2, a1])


def eqnUpperPairs(EIn,n,g):
    """

    :param EIn:
    :param n:
    :param g:
    :return:
    """

    E=EIn
    upperPairsAll=retUpperPair(g,E)
    retValsCis=[]# in the order x2, x1
    retValsTrans=[]#in the order x1, x2

    #######
    #fill cis
    for pairTmp in upperPairsAll:
        x2Tmp,x1Tmp=pairTmp
        try:
            intValTmp=integralQuadrature(g,E,x1Tmp,x2Tmp)
        except ZeroDivisionError:
            continue
        rstTmp = intValTmp - (n + 1 / 2) * mp.pi
        retValsCis.append(rstTmp)
    # # fill trans
    for pairTmp in upperPairsAll:
        x2Tmp, x1Tmp = pairTmp
        try:
            intValTmp = integralQuadrature(g, E, x2Tmp, x1Tmp)
        except ZeroDivisionError:
            continue
        rstTmp = intValTmp - (n + 1 / 2) * mp.pi
        retValsTrans.append(rstTmp)
    ########################################
    retCombined=retValsCis+retValsTrans
    retSorted=sorted(retCombined,key=mpmath.fabs)
    root0=retSorted[0]
    #return complex root
    return root0
    #real
    # return mpmath.fabs(root0)

def computeOneSolutionWithUpperPairs(inData):
    """

    :param inData: [n,g,Eest]
    :return: [n,g,E]
    """
    n, g, Eest = inData
    Eest *= (1 + 0.1j)
    try:
        E = mpmath.findroot(lambda EVal: eqnUpperPairs(EVal, n, g), Eest, solver="muller", maxsteps=100,
                            tol=1e-10)
        return [n, g, E]
    except ValueError as e:
        # print("error message:"+str(e))
        return []





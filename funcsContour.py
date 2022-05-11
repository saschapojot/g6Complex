import numpy as np
import cmath
import mpmath
from mpmath import mp
mp.dps=15

#this script contains functions for contour integral

def retUpperPair(g,E):
    """

    :param g:
    :param E:
    :return: the pair of roots that lie on y>0
    """
    coef = [g, -1j, 0, 0, np.complex128(E)]
    upperPair = []
    rootsAll = np.roots(coef)
    for rt in rootsAll:
        if np.imag(rt) > 0:
            upperPair.append(rt)
    #0th elem  x<0, 1st elem x>0
    upperPair = sorted(upperPair, key=np.real)
    return upperPair

def minusQ(x,g,E):
    return E-1j*x**3+g*x**4

def Qpp(x,g):
    return 6*1j*x-12*g*x**2

def s0s2p(x,g,E):
    return cmath.sqrt(minusQ(x,g,E))+1/48*Qpp(x,g)/cmath.sqrt(minusQ(x,g,E)**3)



def evalIntegral(g,E):
    """

    :param g:
    :param E:
    :return: contour integral
    """
    z2,z1=retUpperPair(g,E)
    r = 0.01
    angle = 0.1 * np.pi
    #coordinates of points on contour
    x1 = np.real(z1)
    y1 = np.imag(z1)
    x2 = np.real(z2)
    y2 = np.imag(z2)

    xA = x1 - r * np.cos(angle)
    xB = x2 + r * np.cos(angle)
    xC = x2
    xD = x1

    yA = y1 + r * np.cos(angle)
    yB = y2 + r * np.sin(angle)
    yD = y1 - r * np.sin(angle)
    yC = y2 - r * np.sin(angle)


    #dx=1e-4
    #number of points on horizontal line
    #M=int(np.abs(x1-x2)/dx)
    # lower line CD
   # CD = 1j * yC + np.linspace(x2, x1, M)
    #zCDVals = [s0s2p(x,g,E) for x in CD]

    #number of points on arc
    N=500

    #right arc
    DArcAThetaVals = np.linspace(angle - np.pi, np.pi - angle, N)
    zDArcAVals = [s0s2p(z1 + r * np.exp(1j * th),g,E) for th in DArcAThetaVals]
    dTheta = DArcAThetaVals[1] - DArcAThetaVals[0]

    #upper line AB
    #AB = 1j * yA + np.linspace(xA, xB, M)
    #zABVals = [s0s2p(x,g,E) for x in AB]

    #arc BC
    BArcCThetaVals = np.linspace(angle, 2 * np.pi - angle, N)
    zBArcCVals = [s0s2p(z2 + r * np.exp(1j * th),g,E) for th in BArcCThetaVals]

    #find jump of DA
    DArcAImag = np.imag(zDArcAVals)
    derivDArcAImag = np.abs(np.diff(DArcAImag)) / dTheta
    inds = np.argsort(derivDArcAImag)[::-1]
    ind1 = inds[0] + 1#the 0th point of jump on arc DA

    #find jump of BC
    BArcCImag = np.imag(zBArcCVals)
    derivBArcCImag = np.abs(np.diff(BArcCImag)) / dTheta
    indsBArcC = np.argsort(derivBArcCImag)[::-1]

    ind2 = indsBArcC[0] + 1#0th point of jump on arc BC

    intCD = mpmath.quad(lambda x: s0s2p(1j * yC + x,g,E), [xC, xD])

    intDAPart1 = mpmath.quad(lambda th: r * 1j * mpmath.exp(1j * th) * s0s2p(z1 + r * mpmath.exp(1j * th),g,E),
                             [DArcAThetaVals[0], DArcAThetaVals[ind1 - 1]])

    # -1
    intDAPart2 = mpmath.quad(lambda th: r * 1j * mpmath.exp(1j * th) * s0s2p(z1 + r * mpmath.exp(1j * th),g,E),
                             [DArcAThetaVals[ind1], DArcAThetaVals[-1]])

    # -1
    intAB = mpmath.quad(lambda x: s0s2p(1j * yA + x,g,E), [xA, xB])

    # -1
    intBCPart1 = mpmath.quad(lambda th: r * 1j * mpmath.exp(1j * th) * s0s2p(z2 + r * mpmath.exp(1j * th),g,E),
                             [BArcCThetaVals[0], BArcCThetaVals[ind2 - 1]])

    intBCPart2 = mpmath.quad(lambda th: r * 1j * mpmath.exp(1j * th) * s0s2p(z2 + r * mpmath.exp(1j * th),g,E),
                             [BArcCThetaVals[ind2], BArcCThetaVals[-1]])
    tmp = intCD + intDAPart1 - intDAPart2 - intAB - intBCPart1 + intBCPart2
    tmp *= 1 / 2
    return tmp


def computeOneSolutionWithUpperPairs(inData):
    """

    :param inData: [n,g,Eest]
    :return: [n,g,E]
    """
    n,g,Eest=inData
    Eest*=1+0j
    try:
        E=mpmath.findroot(lambda EVal: evalIntegral(g,EVal)-(n+1/2)*mpmath.pi,Eest,solver="muller",maxsteps=100,tol=1e-7)
        return [n,g,E]
    except ValueError as e:
        print(str(e))
        return []




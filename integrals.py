import numpy as np
from display import *

initTemp = 100
dT = 50          
temperature = 1200      # temp otoczenia
conducitivity = 25      # przewodność
alfa = 300              # wsp konwekcji [W/m2K]
specificHeat = 700      # cieplo wlasciwe
ro = 7800               # gestosc

def f(x):
    return 5*pow(x,2) + 3*x + 6

def f3D(punkt):
    return 5*pow(punkt[0],2)*pow(punkt[1],2) + 3*punkt[0]*punkt[1] + 6

def calka2D(X, W):
    integral = 0.0
    for i in range(len(X)):
        integral += f(X[i]) * W[i]
    return integral

def getX(punkt):
    return punkt[0]
def getY(punkt):
    return punkt[1]

def matrixSum(mat1, mat2):
    SUM = []
    for i, x in enumerate(mat1):
        row = []
        for j, y in enumerate(mat1):
            row.append( mat1[i][j] + mat2[i][j] )
        SUM.append(row)
    return SUM  

def vectorByTransposition(Arr):
    XX = []
    for x in Arr:
        row = []
        for y in Arr:
            row.append(x*y)
        XX.append(row)

    return XX

def matrixByNumber(matrix, det):
    H = []
    
    for row in matrix:
        ROW = []
        for number in row:
            ROW.append(number * det)
        H.append(ROW)
    
    return H

def punktydo3D(X):
    lp = len(X)
    punkty = []

    for i in range(0,len(X)):  
        for j in  range(0,len(X)):
            punkty.append( [X[i], X[j]] )
    punkty.sort()
    punkty.sort(key=getY)

    for i in range(lp, len(punkty)-lp+1):
        if i%lp == 0:
            punkty[i], punkty[i+lp-1] = punkty[i+lp-1], punkty[i]
        

    # for i, pkt in enumerate(punkty):
    #     if i<len(punkty)-1:
    #         if punkty[i][1] == punkty[i+1][1] and punkty[i][0] < punkty[i+1][0]:
    #             punkty[i], punkty[i+1] = punkty[i+1], punkty[i]
    return punkty

def calka_3D(punkty, W):        
    integral = 0.0
    Wn = []

    # lista kolejnych iloczynow wag (w1*w2...)
    for i in range(0,len(W)):  
        for j in  range(0,len(W)):
            Wn.append(W[i] * W[j])

    # obliczenie calki mnozac iloczyn wag przez wartosc funkcji
    for i, pkt in enumerate(punkty):
        integral += Wn[i] * f3D(pkt)

    return integral
    
def calkowanieTest():
    # print("1d 1 punktowy")c
    # print(calka2D(N1x, N1w))
    # print("1d 2 punktowy")
    # print(calka2D(N2x, N2w))
    # print("1d 3 punktowy")
    # print(calka2D(N3x, N3w))

    #pkt = punktydo3D(N1x)
    # print(pkt)
    # print("2d 2 punktowy")
    # print(calka_3D(pkt, N1w))
    # print("2d 3 punktowy")
    # print(calka_3D(pkt1, N2w))
    #xi = xiArray(pkt)

    #eta = etaArray(pkt)
    #print(eta)
    return 0


# def xiFunc(point):
#     waga = (8/9 if point[1] == 0.0 else 5/9)

#     return [(-1.0 / 4.0) * (1 - point[1]) * waga, (1.0 / 4.0) * (1 - point[1]) * waga, (1.0 / 4.0) * (1 + point[1]) * waga,  (-1.0 / 4.0) * (1 + point[1]) * waga]


# def etaFunc(point):
    waga = (8/9 if point[1] == 0.0 else 5/9)
    return [(-1.0 / 4.0) * (1 - point[0]) * waga, (-1.0 / 4.0) * (1 + point[0])* waga, (1.0 / 4.0) * (1 + point[0])* waga,  (1.0 / 4.0) * (1 - point[0]) * waga]

def xiFunc(point):
    return [(-1.0 / 4.0) * (1 - point[1]), (1.0 / 4.0) * (1 - point[1]), (1.0 / 4.0) * (1 + point[1]),  (-1.0 / 4.0) * (1 + point[1])]


def etaFunc(point):
    return [(-1.0 / 4.0) * (1 - point[0]), (-1.0 / 4.0) * (1 + point[0]), (1.0 / 4.0) * (1 + point[0]),  (1.0 / 4.0) * (1 - point[0])]


def shapeFunc(point):
    '''
    Parameter:
        point (list): coordinates of which are given into formula
    Formulas for shape functions in 4-node element
    '''
    return [ 0.25 * (1 - point[0]) *(1 - point[1]), 0.25 * (1 + point[0]) * (1 - point[1]), 0.25 * (1 + point[0])*(1 + point[1]), 0.25 * (1 - point[0]) * (1 + point[1]) ]


def xiArray(points):
    elem = []
    for i,point in enumerate(points):
        elem.append(xiFunc(point))
    return elem


def etaArray(points):
    elem = []
    for i,point in enumerate(points):
        elem.append(etaFunc(point))
    return elem


def calculateJacobian(element, n, xi, eta):
    '''
    Description:
    Calculates Jacobian matrix, given xi and eta matrices,
    n is the following row index in xi/eta array
    '''
   
    J=[] 
    dYdEta = 0.0
    dYdXi = 0.0

    dXdEta = 0.0
    dXdXi = 0.0
    
    xList = []
    yList = []
    for ID in element.IDcoords:
        xList.append(ID[0])
        yList.append(ID[1]) 
        
    
    for i,x in enumerate(xList):
        dXdXi += x * xi[n][i]
        dXdEta += x * eta[n][i]
        #
        # ((x,xi[0][i], dXdXi))

    for i,y in enumerate(yList):
        dYdEta += y * eta[n][i]
        dYdXi += y * xi[n][i]
        
    J.append(dYdEta)
    J.append(-dYdXi)
    J.append(-dXdEta)
    J.append(dXdXi)
        
    return J


def jakobianInv(J):
    #print("1/det:")
    d = det(J)
    #print(1/d)
    

    for i,index in enumerate(J):
        J[i] *= 1/d
   # print(det(J))
    return J


def det(J):
    return J[0]*J[3] - J[1]*J[2]


def arrH(Arr, n):
    '''
    Returns n-th row in array multiplied by its transposition
    '''
    Row = Arr[n]
    return vectorByTransposition(Row)


def matrixH(SUM, conducitivity, det):
    '''
    Parameters:
        SUM (list): sum matrix of dX and dY, arrays
        conducitivity (double): conductivity of material
        det (double): determinant of J matrix
    
    Returns:
        2D list, H matrix, that shows how much heat given element can transport
    '''
    H = []
    for i, x in enumerate(SUM):
        row = []
        for j, y in enumerate(SUM):
            row.append( SUM[i][j] * det * conducitivity )
        H.append(row)
    return H  

def matrixC(NNt, Heat, ro, det):
    '''
    Parameters:
        NNt (list): array made by multiplication of shape func vector and traspositon of that vector
        Heat (double):  specific heat of material
        ro (double): density of material
        det (double): determinant of J matrix
    
    Returns:
        2D list, C matrix, that shows how much heat given element can accumulate
    '''
    C = []
    for i, x in enumerate(NNt):
        row = []
        for j, y in enumerate(NNt):
            row.append( NNt[i][j] * det * ro * Heat *dT )
        C.append(row)
    return C





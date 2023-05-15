from elem4_2D import *

class SOE:
    def __init__(self, grid, elem) :
        self.grid = grid
        self.globalH = setGlobalH(self.grid)
        self.globalH_HBC = setGlobalH_HBC(self.grid)
        self.globalC = setGlobalC(self.grid)
        self.globalP = elem.vecP
        self.H_C_dT = self.calculateH_C_dT()
        self.el = elem
    
    def calculateH_C_dT(self):
        '''
        Returns:
        Sum of C/dT and global[H+HBC] matrices
        '''
        C = self.globalC
        C_dT = []

        for i,x in enumerate(self.globalC):
            row = []
            for j,y in enumerate(self.globalC):
                C[i][j] = C[i][j]/dT
                row.append( C[i][j]/dT)
            C_dT.append(row)
        result = self.globalH_HBC
        #displayArr(C, "C/dT")
        return matrixSum(result, C_dT)

    def calculateP_C_dT(self, I):
        '''
        Description:
        Calculates vector {P}+{[C]/dT}*{T0} and PRINTS min max T values in time
        after solving equation of [H]+[C]/dT and above vector
        '''
        vector = np.zeros(self.grid.nN)
        
        for i in range(len(self.globalP)):
            vector[i] = self.globalP[i]
        
        globalC = np.mat(self.globalC)

        for i in range(self.grid.nN):
            for j in range(self.grid.nN):
                # {P} = {P}+{[C]/dT}*{T0}
                vector[i] +=  (globalC[i,j] * self.grid.nodes[j].t0 / dT ) 
        #print(I, vector)
        from scipy import linalg
        tempVec = linalg.solve(self.H_C_dT, vector)

        for x in range(len(tempVec)):
            # updating temperatures in nodes
            self.grid.nodes[x].t0 = tempVec[x]
            #print(t, self.grid.nodes[x].t0 )
        #displayArr(self.H_C_dT,"XD")
        vec = tempVec
        print(I*dT,"\t", np.amin(tempVec),"\t", np.amax(tempVec))
        
        

    def displayResults(self):
        print("{x} points of integration system".format(x = len(self.el.W)) )
        print("Time[s]\t MinTemp[°C]\t\t MaxTemp[°C]")
        for i in range(1,21):
            self.calculateP_C_dT(i)

def setGlobalH(grid):
    '''
    Description:
    Local H matrices in elements are aggregated into one global matrix
    '''
    globalH = []

    for i in range(grid.nH * grid.nB):
        row = []
        for j in range(grid.nH * grid.nB):
            row.append(0)
        globalH.append(row)

    for i,element in enumerate(grid.elements):
        for x,ID1 in enumerate(element.ID):
            for y,ID2 in enumerate(element.ID):
                globalH[ID1-1][ID2-1] += element.Hsum[x][y] 
    
    return globalH


def setGlobalH_HBC(grid):
    '''
    Description:
    Local H + HBC matrices in elements are aggregated into one global matrix
    '''
    globalHBC = []

    for i in range(grid.nH * grid.nB):
        row = []
        for j in range(grid.nH * grid.nB):
            row.append(0)
        globalHBC.append(row)

    for i, element in enumerate(grid.elements):
        for x,ID1 in enumerate(element.ID):
            for y,ID2 in enumerate(element.ID):
                globalHBC[ID1-1][ID2-1] += element.H_Hbc[x][y] 
    
    return globalHBC

def setGlobalC(grid):
    '''
    Description:
    Local C matrices in elements are aggregated into one global matrix
    '''
    globalC = []

    for i in range(grid.nH * grid.nB):
        row = []
        for j in range(grid.nH * grid.nB):
            row.append(0)
        globalC.append(row)

    for i, element in enumerate(grid.elements):
        for x,ID1 in enumerate(element.ID):
            for y,ID2 in enumerate(element.ID):
                globalC[ID1-1][ID2-1] += element.Csum[x][y] 
    
    return globalC


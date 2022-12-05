from integrals import *
from display import displayArr


class elem4_2D:
    '''
        Description:
            Object represents universal 4-node 2-dimensional element

        Args:
            - grid (GRID): grid of finite elements
            - points (list): universal element points
            - xi (list): array of dN/dXi functions
            - eta (list): array of dN/dEta functions
            - Nx (list): integration scheme nodes
            - Nw (list): integration scheme coefficients
            - W (list): integration scheme weights
            - surfacePoints (list): universal element surface points
            - N (list): shape function array for surfacePoints
            - NC (list): shape function array for point
            - weightsH (list): multiplicated corresponding weights (used to calc matrixH and matrixC)
            - jakobian (list): jakobian array 
            - walls (list): list of surface walls (3 surface points = 1 wall) 
            - vecP (list): all P vectors aggregated into one vector  


    '''
    def __init__(self, grid, Nx, Nw, W):
        '''
        Args:
            - grid (GRID): finite elements grid 
            - Nx (list): integration scheme nodes
            - Nw (list): integration scheme coefficients
            - W (list): integration scheme weights
        '''
        
        self.grid = grid
        self.points = punktydo3D(Nx)
        self.xi = xiArray(self.points)
        self.eta = etaArray(self.points)
    
        self.Nx = Nx
        self.Nw = Nw
        self.W =  W
        
        self.surfacePoints = self.SurfacePoints()

        self.N = self.N()
        self.NC = self.NC()
        self.weightsH = self.WeightsH()
        
        self.setElementsMatrixH()
        self.jakobian = self.JakobianElem()

        self.walls = self.Walls()
        self.setElementsHbc()

        self.setElementsVectorP()
        self.vecP = self.AggregatedP()
        self.setElementsMatrixC()
        self.TestElemMatrixH()

    def JakobianElem(self):
        '''
        Description:
        Returns Jacobian containing derivatives of shape functions
        '''
        J=[] 
        dYdEta = 0.0
        dYdXi = 0.0

        dXdEta = 0.0
        dXdXi = 0.0
      
        for i,point in enumerate(self.points):
            for j in range(4):
                dXdXi += point[0] * self.xi[0][j]
                dXdEta += point[0] * self.eta[0][j]

                dYdEta += point[1] * self.eta[0][j]
                dYdXi += point[1] * self.xi[0][j]

        J.append(dYdEta)
        J.append(-dYdXi)
        J.append(-dXdEta)
        J.append(dXdXi)
            
        return J


    def TestElemMatrixH(self):
        Hsum = []
        
        # [[0, 0, 0, 0],
        #    ...
        #  [0, 0, 0, 0]]    
        J = self.jakobian
        xMultiply = J[0]
        yMultiply = J[3]

        xArray = matrixByNumber(self.xi, xMultiply)
        yArray = matrixByNumber(self.eta, yMultiply)
    
        for j in range(4):
            row = [0 for i in range(4)]
            Hsum.append(row)
            

        for j, point in enumerate(self.points):
            HX = arrH(xArray, j)
            HY = arrH(yArray, j)
            S = matrixSum(HX,HY)
            H = matrixH(S, conducitivity, 6400)
            
           
            Hsum = matrixSum(Hsum, H)
        
    def SurfacePoints(self):
        '''
        Description:
            Returns points placed on surface of element 
        '''
        sPoints = []
    
        for p in self.Nx:
            for n in range(-1,2,2):
                sPoints.append([p, n])
                sPoints.append([n, p])


        sPoints = sorted(sPoints, key=lambda x: x[0])
        #sPoints = sorted(sPoints, key=lambda x: x[1])
        return sPoints

    def N(self):
        '''
        Description:
            Returns shape functions list of points on surface
        '''
        N = []
        for point in self.surfacePoints:
            N.append(shapeFunc(point))
        
        return N

    def NC(self):
        '''
        Description:
            Returns shape functions list of all points
        '''
        NC = []
        for point in self.points:
            NC.append(shapeFunc(point))
        #displayArr(NC,"XDD")
        return NC

    def testC(self):
        allC = []
        for j in self.points:
                row = [0 for i in self.points]
                allC.append(row)

        for row in self.NC:
            NNt = []
            for x in row:
                r = []
                for y in row:
                    r.append(x * y * ro * specificHeat * 0.00015625)
                NNt.append(r)
            allC = matrixSum(allC, NNt)
            displayArr(allC,"C")

    def Walls(self):
        '''
        Description:
            Returns 'coords' of universal element surface walls points, that are vertices of thet element
        '''
        walls = []
        minCoord = min(min(self.surfacePoints))
        maxCoord = max(max(self.surfacePoints))

        walls.append([ p for p in self.surfacePoints if p[0] == minCoord ]) 
        walls.append([ p for p in self.surfacePoints if p[1] == minCoord ])
        walls.append([ p for p in self.surfacePoints if p[0] == maxCoord ])
        walls.append([ p for p in self.surfacePoints if p[1] == maxCoord ]) 

        SORTED = []
        #print("walls")
        for wall in walls:
            wall = sorted(wall,key=lambda l:l[1])
            SORTED.append(wall)
            #print(wall)
        return SORTED
    
    def WeightsH(self):
        '''
        Description:
            Returns weights list used in two-dimensional operations
        '''
        wagi = []
        for x in self.W:
            for y in self.W:
                wagi.append(x*y)
        return wagi

    def setElementsMatrixH(self): 
        '''
        Description:
            Sets H matrices for all elements in grid, H matrix is representation Fourier-Kirchoff equation
        '''
        for i,el in enumerate(self.grid.elements):
            Hsum = []
          
            for j in el.ID:
                row = [0 for i in el.ID]
                Hsum.append(row)

            for j, point in enumerate(self.points):
                el.jakobian = calculateJacobian(el,j, self.xi, self.eta)
                
                J = el.jakobian

                Jinv = jakobianInv(J)
            
                xMultiply = el.jakobian[0]
                yMultiply = el.jakobian[3]

                xArray = matrixByNumber(self.xi, xMultiply  )
                yArray = matrixByNumber(self.eta, yMultiply )

                HX = arrH(xArray, j)   
                HY = arrH(yArray, j)
                S = matrixSum(HX,HY)
                el.H =matrixH(S, conducitivity, 1/det(J) * self.weightsH[j])
                el.P = []
                el.Hlist.append(el.H)
                Hsum = matrixSum(Hsum, el.H)
            el.Hsum = Hsum

    def setElementsHbc(self):
        '''
        Description:
            Sets Hbc matrices for all elements in grid, Hbc represents boundary condition integral in [H] matrix equation
        '''
        for i,el in enumerate(self.grid.elements):
            SUM = []
            for j in el.ID:
                row = [0 for i in el.ID]
                SUM.append(row)
            
            for j, wall in enumerate(el.wallsBC):
                el.jakobian = calculateJacobian(el,j, self.xi, self.eta)
                J = el.jakobian
                
                if wall == 1:                                                       # if element wall on boundary
                    elemHBC = calculateHbc(self.walls, self.grid.deltaX/2, self.W)  # calculate all Hbc matrices in universal element
                    el.HbcWalls.append(elemHBC[j])                                  # append right matrix to grid element all Hbc list   

            el.HbcSum = el.sumHbc()                                                 # sums all Hbc matrices in element
            SUM = matrixSum(el.Hsum, el.HbcSum)                                     # sums H and Hbc sum matrices
            el.H_Hbc = SUM                                                          # assigns result to element

    def setElementsMatrixC(self):
        '''
        Description:
            - Sets C matrices for elements in grid
            - Calls matrixC func
        '''
        for i,el in enumerate(self.grid.elements):
            Csum = []
            
            # [[0, 0, 0, 0],
            #    ...
            #  [0, 0, 0, 0]]    
            for j in el.ID:
                row = [0 for i in el.ID]
                Csum.append(row)
         
            for j, node in enumerate(self.points):
                J = el.jakobian
                NNt = []
                for x in self.NC[j]:
                    row = []
                    for y in self.NC[j]:
                        row.append(x*y)
                    NNt.append(row)    
                C = matrixC(NNt, specificHeat, ro, det(J) * self.weightsH[j]  )
                Csum = matrixSum(Csum, C)
            el.Csum = Csum

    def setElementsVectorP(self):
        '''
        Description:
            Calls calculateVectorP for every element in grid
        '''
        for el in self.grid.elements:
           el.vectorP = calculateVectorP(el, self.walls, self.grid.deltaX/2, self.W)

    def AggregatedP(self):
        '''
        Description:
            Calls aggregateVectorP for every element in grid,
            then returns result of these aggregations 
        '''
        aggregatedP= [0 for i in range(self.grid.nH * self.grid.nB)]
        for el in self.grid.elements:
            aggregateVectorP( el, aggregatedP)
        return aggregatedP


def calculateHbc(walls, det, W):
    '''
    Description:
        - Creates Hbc matrix given element walls, det[J] and weights
        - Returns all of Hbc matrices in universal element 
    '''
    allHBC = []
    for i,wall in enumerate(walls):
        
        shapes = []
        for x in range(len(wall)):       
           # c
            shapes.append( shapeFunc(wall[x]))
        
        SUM = [[0,0,0,0],
               [0,0,0,0],
               [0,0,0,0],
               [0,0,0,0]]
        #print()
        for i,shape in enumerate(shapes):
            SUM = matrixSum(SUM, matrixByNumber(vectorByTransposition(shape), W[i]))
    
        HBC = matrixByNumber(SUM, det * alfa  )
        allHBC.append(HBC)
    return allHBC


def calculateVectorP(element, walls, det, W):
    '''
    Description:
        Returns Vector P for given element
            1) calculated for points on surface of element
            2) sums all vectors on every BC wall
    '''
    allP = []
    #print(len(walls))
    for i,wall in enumerate(walls):
        if element.wallsBC[i] == 1:                     # if wall on boundary
            shapes = []
            for x in range(len(wall)):       
                shapes.append( shapeFunc(wall[x]))      # append shape function for every point in wall

            SUM = [ 0 for d in range(len(shapes[0]))]
            for i,shape in enumerate(shapes): 
                for j,func in enumerate(shape):
                    SUM[j] += func * temperature * det * alfa * W[i] # formula
            allP.append(SUM)
    #print(allP)
    X= [ 0 for i in range(4) ]

    # sums P vectors from every wall
    for P in allP:
        for i, number in enumerate(P):
            X[i] += number
    return X


def aggregateVectorP( element, aggP):
    '''
    Description:
    Changes indices of aggP depending on IDs of element
    '''
    for i,ID in enumerate(element.ID):
        aggP[ID-1] += element.vectorP[i]
    
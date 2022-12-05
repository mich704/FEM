
from integrals import *

class node():
    '''
        Args:
            - BC (int): node boundary condition (1- boundary, 0- 'inside' grid)
            - t0 (double): node initial temperature
    '''
    def __init__(self, x, y, id):
        '''
            Args:
                - x (double): horizontal coordinate
                - y (double): vertical coordinate
                - id (int): identificator
        '''
        self.id = id
        self.x = x
        self.y = y
    
        self.BC = 0.0
        self.t0 = initTemp 

    def __str__(self):
        return str(str(self.id)+": "+str((self.x, self.y)))

class element():
    '''
        Description: 
            Represents single element in FEM grid
        
        Args:
            - ID (list): element IDs in numbering order
            - IDcoords (list): element IDs coordinates in grid layout
            - nodes (list): list of nodes objects within element
            - jakobian (list): Jacobi matrix, 
            - Hlist (list): H (Fourier-Kirchoff equations) matrices list
            - Hsum (list): sum of matrices within Hlist
            - wallsBC (list): element walls boundary conditions
            - HbcWalls (list): Hbc matrices list
            - HbcSum (list): sum of matrices within HbcWalls
            - H_Hbc (list): matrix sum of HbcSum and Hsum 
            - vectorP (list): element's vector P
            - Csum (list): sum of C matrices within element

    '''

    def __init__(self, id1, grid):
        '''
            Args:
                - id1 (int): first id of element
                - grid (GRID): grid object
            
        '''
        self.ID = []

        self.id1 = id1
        self.id2 = id1 + grid.nH
        self.id3 = self.id2 + 1
        self.id4 = self.id1 + 1
        self.grid = grid

        self.ID.append(self.id1)
        self.ID.append(self.id2)
        self.ID.append(self.id3)
        self.ID.append(self.id4)

        self.IDcoords = []
        self.nodes = []
        self.jakobian = []
    
        self.Hlist = []
        self.Hsum = []
    
        self.wallsBC = []
        self.HbcWalls = []
        self.HbcSum = []
        self.H_Hbc = []

        self.vectorP = []
        self.Csum = []
        
        self.setIDCoords()
        self.setNodes()
            
    def setNodes(self):
        nodex = []
        for node in self.grid.nodes:
            if(node.id in self.ID):
                nodex.append(node)
        self.nodes.append(nodex[0])
        self.nodes.append(nodex[2])
        self.nodes.append(nodex[3])
        self.nodes.append(nodex[1])
    
    def setIDCoords(self):
          for ID in self.ID:
            self.IDcoords.append(self.grid.nodesCoords[ID-1])

    def sumHbc(self):
        SUM = []
        for j in self.ID:
            row = [0.0 for i in self.ID]
            SUM.append(row)

        for matrix in self.HbcWalls:
            SUM = matrixSum(SUM, matrix)
        
        return SUM
    
    def __str__(self):
        return str(self.ID)


class GRID():
    def __init__(self, H, B, nH, nB):
        '''
            Args:
                - H (double): vertical dimension 
                - B (double): horizontal dimension
                - nH (int): vertical number of nodes
                - nB (int): horizontal number of nodes  
            
        '''
        self.H = H 
        self.B = B
        self.nH = nH
        self.nB = nB

        self.nN = self.nH * self.nB
        self.xE = self.nB - 1 
        self.yE = self.nH -1 
        self.nE = ( self.xE ) * ( self.yE )
        
        self.deltaY = self.H / ( self.nH -1 )
        self.deltaX = self.B / ( self.nB -1 )

        self.nodes = self.NodesBoard()
        self.elements = []
        self.nodesCoords = self.NodesCoords()
        
        self.setBCFlags()
        self.setElements()
        self.setWallsBC()
        
    def FirstIDs(self):
        '''
            Returns elements first IDs list
        '''

        firstIDs = []
        mod = 0

        for i in range(1, self.nE + 1):
            if i != 1:
                if i % self.yE == 1:
                    mod+=1
                firstIDs.append(i+mod)
            else:
                firstIDs.append(i)

        return firstIDs

    def setBCFlags(self):
        '''
            Sets nodes boundary condition flags
        '''
        for node in self.nodes:
            if(node.x * node.y == 0.0 or node.x == self.B or node.y == self.H):
                node.BC = 1
            else:
                node.BC = 0
        
    def NodesBoard(self):
        '''
            Returns nodes coordinates on grid scale
        '''                          
        XX = 0.0
        YY = 0.0
        
        nodes_coords = []
    
        for el in range(1,self.nN+1):
            n = node(XX,YY,el)
            nodes_coords.append(n)
            YY += self.deltaY
           
            if YY>self.H:
                YY = 0.0
                XX += self.deltaX

        return nodes_coords

    def setElements(self):
        '''
            Set element class objects list
        '''
        firstIDs = self.FirstIDs()

        for i in firstIDs:
            e = element(i, self)
            self.elements.append(e)

    def showElements(self):
        '''
            Prints elements
        '''
        i = 1
        for e in self.elements:
            print("element "+str(i)+"\t"+str(e))
            i += 1

    def showNodes(self):
        '''
            Prints nodes
        '''
        i = 1
        coords = []
        for n in self.nodes:
            print("node "+str(n))
            i += 1

    def NodesCoords(self):
        '''
            Returns nodes coordinates 
        '''
        coords = []
        for node in self.nodes:
            coords.append([node.x, node.y])
        return coords

    def setWallsBC(self):
        '''
            Sets element walls boundary codition flags 
        '''
        BCS=[]
        for element in self.elements:
            element.wallsBC.append(1 if (element.nodes[0].BC == 1 and element.nodes[3].BC == 1) else 0)
            element.wallsBC.append(1 if (element.nodes[0].BC == 1 and element.nodes[1].BC == 1) else 0)
            element.wallsBC.append(1 if (element.nodes[1].BC == 1 and element.nodes[2].BC == 1) else 0)
            element.wallsBC.append(1 if (element.nodes[2].BC == 1 and element.nodes[3].BC == 1) else 0)
            
    # def show_IDcoords(self):
    #     for el in self.elements:
    #         for ID in el.ID:
    #             el.pointID(ID)
    
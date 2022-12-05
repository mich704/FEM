from integrals import *
import matplotlib.pyplot as plt
import numpy as np


def drawPlot(points):
    x= []
    y= []
    #pkt1 = punktydo3D(N2x)
    
    for p in points:
        x.append(p[0])
        y.append(p[1])

    x_coordinates = x
    y_coordinates = y
    plt.scatter(x_coordinates, y_coordinates)
    plt.xlabel('x - axis')
    # naming the y axis
    plt.ylabel('y - axis')
    
    # giving a title to my graph
    plt.title('Grid')
    for i,p in enumerate(points):
        plt.annotate(str(i+1), p)
    plt.savefig('grid.png')
    plt.show()


def displayArr(Arr, name):
    print(name)
    data = np.asarray(Arr)

    #2D array
    if len(data.shape) ==2:
        for row in Arr:
            ROW = []
            for number in row:
                if number == 0.0:
                    ROW.append(0)
                else:
                    ROW.append(round(number,4))
            print(ROW)
       
    #1D array
    if len(data.shape) ==1:
        print(*Arr)
    print()        


def displayElementsMatrices(grid):
    for i,el in enumerate(grid.elements):
        print("element {}".format( i+1))

        for HBC in el.HbcWalls:
            displayArr(HBC, "HBC")
        displayArr(el.HbcSum, "HBC SUM")

        #print()
        for j, node in enumerate(el.nodes):
    
            print()
            
            print("npc: ", j+1)
            #displayArr(xArray, "dN/dX:")

            #displayArr(yArray, "dN/dY:")
            #print(el.jakobian)
            #print(det(el.jakobian))
            #displayArr(node.HX, "HX:")
            #displayArr(node.HY, "HY:")
            #displayArr(matrixSum(node.HX, node.HY), "Hsum {}:".format(j+1))
            
        
        #displayArr(el.Hsum, "Hsum:")            
        print("-------------------------")

def displayElementsHbc(grid):
      for i,el in enumerate(grid.elements):
        print("element {}".format( i+1))  
        for j,wall in enumerate(el.HbcWalls):
            displayArr(wall, "Hbc {}".format(j+1))
        displayArr(el.HbcSum, "Hbc Sum")
        print("-----------------------------------")

def displayElementsH(grid):
      for i,el in enumerate(grid.elements):
        print("element {}".format( i+1))  
        for j,wall in enumerate(el.Hlist):
            displayArr(wall, "H {}".format(j+1))
        displayArr(el.Hsum, "H Sum")
        print("-----------------------------------")

def displayElementsP(grid):
      for i,el in enumerate(grid.elements):
        print("element {}".format( i+1))  
        displayArr(el.vectorP, "vector P")
        print("-----------------------------------")

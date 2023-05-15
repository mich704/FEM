from math import *
from GRID import *
from SOE import *
from display import *
from integrals import *

def main():
    N1x = [-1/sqrt(3), 1/sqrt(3)]
    N1w = [-1,1]
    W1 = [1,1]

    N2x = [-sqrt(3/5), sqrt(3/5), 0]
    N2w = [5/9, 5/9, 8/9]
    W2 = [5/9, 8/9, 5/9]
   

    N3x = [-0.861136, 0.861136, -0.339981, 0.339981 ]
    N3w = [0.347855, 0.347855, 0.652145, 0.652145]
    W3 =  [0.347855, 0.652145,  0.652145 , 0.347855]

    # grid params   
    H = 0.1
    B = 0.1
    nH = 4
    nB = 5
    g = GRID(H, B, nH, nB)
    # print(g.nodesCoords)

    testx = [0.0, 0.025]
    testw =[1.0, 1.0]
    testowe = punktydo3D(testx)
    
    

    print("Please choose integration schema (2, 3 or 4 points): ")
    schemat = int(input())

    if schemat == 2:
        elem2 = elem4_2D(g, N1x, N1w, W1)
        soe = SOE(g, elem2)
    
        soe.displayResults()

    elif schemat == 3:
        elem3 = elem4_2D(g, N2x, N2w, W2)
        soe = SOE(g, elem3)
        soe.displayResults()

    elif schemat == 4:
        elem3 = elem4_2D(g, N3x, N3w, W3)
        soe = SOE(g, elem3)
        soe.displayResults()

    else:
        print("nie ma takiego schematu")


    print()
    #displayElementsP(g)
    #displayArr(elem.eta,"dN/dEta: ")
    #displayArr(elem.xi,"dN/dXi: ")
    
    #displayElementsMatrices(g)
    #displayElementsHbc(g)
    #displayElementsH(g)
   
    print()
    print()

    # # for el in g.elements:
    # #     displayArr(el.jakobian, "JAKOBIAN")
    # #     #print(1/det(el.jakobian))
    # # print("deltaX", "\tdeltaY")
    # # print(g.deltaX,"\t", g.deltaY)
    # # displayArr(soe.globalH, "global H: ") 
    # # displayArr(soe.globalH_HBC, "global H + HBC: ")
    # # displayArr(soe.globalC, "global C: ")
    # # displayArr(soe.globalP, "vecP")
    # # displayArr(soe.H_C_dT, "[H] = [H]+[C]/dT ")
    # soe.displayResults()
    # #drawPlot(g.nodesCoords)

if  __name__ == '__main__':
    main() 

import numpy as np
import scipy
import scipy.linalg
import sys
from collections import Counter
# %%

""" construc huckel matrix for linear n-polyene """
def polyene_huckel_matrix (n,shape_polyene):
    matrix_polyene = np.zeros(shape=(n,n))
    
    prev = 0
    for i in range(n+1): #n+1 because index starts from 0
        if prev != 0 and i != 0 and i != prev:
            matrix_polyene[i - 1, prev - 1] = -1
            matrix_polyene[prev - 1, i - 1] = -1
        prev = i
    #for cyclic, further modify the 2 corners of the matrix
    if shape_polyene == "cyclic":
        matrix_polyene[0, n-1] = -1
        matrix_polyene[n-1, 0] = -1

    return matrix_polyene

""" hard-coded the matrix for platonic solids """
def platonic_huckel_matrix(shape_platonic):
    if shape_platonic == "tetrahedron":
        matrix_platonic = [[ 0, -1, -1, -1],
                           [-1,  0, -1, -1],
                           [-1, -1,  0, -1],
                           [-1, -1, -1,  0]]

    elif shape_platonic == "cube":
        matrix_platonic = [[ 0, -1,  0, -1,  0,  0,  0, -1],
                           [-1,  0, -1,  0,  0,  0, -1,  0],
                           [ 0, -1,  0, -1,  0, -1,  0,  0],
                           [-1,  0, -1,  0, -1,  0,  0,  0],
                           [ 0,  0,  0, -1,  0, -1,  0, -1],
                           [ 0,  0, -1,  0, -1,  0, -1,  0],
                           [ 0, -1,  0,  0,  0, -1,  0, -1],
                           [-1,  0,  0,  0, -1,  0, -1,  0]]

    elif shape_platonic == "dodecahedron":
        matrix_platonic = [[ 0, -1,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                           [-1,  0, -1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                           [ 0, -1,  0, -1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                           [ 0,  0, -1,  0, -1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                           [-1,  0,  0, -1,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                           [-1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0],
                           [ 0, -1,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0],
                           [ 0,  0, -1,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0, -1,  0,  0,  0,  0,  0],
                           [ 0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0],
                           [ 0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0],
                           [ 0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0],
                           [ 0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0],
                           [ 0,  0,  0,  0,  0, -1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0],
                           [ 0,  0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  -1, 0],
                           [ 0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1],
                           [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0, -1,  0,  0, -1],
                           [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0, -1,  0, -1,  0,  0],
                           [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0, -1,  0,  -1, 0],
                           [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0, -1,  0, -1],
                           [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  -1, 0]]

    return matrix_platonic

""" get Hückel pi-energies """
def get_evals (matrix):
    evals, evects = np.linalg.eig(matrix)
    energies = sorted(round(num,3) for num in evals) #overcome round-off error 
    
    return energies


""" get degeneracy of the energies """
def degeneracy ():
    energies = get_evals(matrix)
    print("The energies of the molecule are:\n" + str(energies)) 
    print("------------------------------")
    _count = Counter()
    _count.update(energies)
    for e in _count.keys():
        print("Denegeracy of energy %s: %d" % (e, _count[e])) 


#get information about the molecule and execute the code 
flag = False
while not flag:
    molecule_type = input("Type of the molecule is 1-polyene; 2-sp2 Platonic solid (enter the number):\n")
    if molecule_type in ["1", "2"]:
        flag = True
    else:
         print("Invalid input. Please try again")


if molecule_type == "1":
    flag = False
    while not flag:
        shape_polyene = input("linear or cyclic:\n")
        n = int(input("number of atoms\n")) 
        if shape_polyene in ["linear", "cyclic"] and n >0:
            flag = True
        else: 
            print ("This is not a valid polyene molecule") 
        
    matrix = polyene_huckel_matrix(n, shape_polyene)  
    get_evals(matrix)
    degeneracy()



elif molecule_type == "2":
    flag = False
    while not flag: 
        shape_platonic = input("tetrahedron or cube or dodecahedron:\n") 
        if shape_platonic in ["tetrahedron", "cube", "dodecahedron"]:
            flag = True
        else:
            print ("Sorry, we do not have information for now. Please enter one of the three platonic solids given.")

    matrix = platonic_huckel_matrix(shape_platonic)
    get_evals(matrix)
    degeneracy()

else:
    print("Sorry, we do not have information for now.")
#%%




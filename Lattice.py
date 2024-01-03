import numpy as np
import matplotlib.pyplot as plt
import copy
import Percolation_functions as func

class Lattice:
    def __init__(self,size: int, p):
        self.size = size
        self.p = p
        self.lat = np.array(np.random.binomial(1, p, (size,size)))
    
    #-------- getters -------
    def __getlat__(self):
        return self.lat
        
    def __getsize__(self):
        return self.size
    
    def __getp__(self):
        return self.p
    # ------ setters --------
    def __setLat__(self, lattice):
        self.lat = lattice
        self.size = lattice.shape[0]

    #-------- methods -------
            
    def burningMethod(self):
        '''Burn that lattice!
        output: 
            burning_forest: NDarray - lattice with increasing site indexes
            isEndreached(...): Bool - duh
            t-1: Int - lenght of shortest path to end row (in case of reached end)'''

        #initial variables
        burning_forest = copy.deepcopy(self.lat)
        size = copy.deepcopy(self.size)
        t=2
        ref_burning_forest = copy.deepcopy(burning_forest)
        #step 1) burning top row
        for i in range(self.size):
            if burning_forest[0][i] == 1:
                burning_forest[0][i] = t

        #while bottomm isnt reached and there is still fire (ref and burning arent equal)
        while not func.isEndReached(burning_forest, size, t) and (not np.array_equal(burning_forest, ref_burning_forest)):
            #actualize reference
            ref_burning_forest = copy.deepcopy(burning_forest)
            #burn every next cell
            for i in range(self.size):
                for j in range(self.size):
                    if burning_forest[i][j] == t:
                        if burning_forest[func.NI(size, i)][j] == 1:  
                            burning_forest[func.NI(size, i)][j] = t + 1
                        if burning_forest[func.PI(i)][j] ==1:
                            burning_forest[func.PI(i)][j] = t+1 
                        if burning_forest[i][func.NI(size,j)]==1:
                            burning_forest[i][func.NI(size,j)]= t+1
                        if burning_forest[i][func.PI(j)]==1: 
                            burning_forest[i][func.PI(j)]= t+1
            #actualize t           
            t = t+1
        return burning_forest, func.isEndReached(burning_forest, size, t), t-1
    
    #Hoshen-Kopelman algorithm
    def HoshenKopelman(self):
        '''Method returning clustered lattice and their sizes based on HK algorythm
        outputs:   
            N_ij: 2D NDarray - clustered lattice
            M_k: 1D NDarray - sizes of clusters primary indexed (with zeros and negative values) '''
        size = self.size
        k = 2
        sizeMK = int(size*size/2)
        M_k = np.zeros(sizeMK, dtype=int)
        N_ij = copy.deepcopy(self.lat)
            
        #labeling clusters
        for i in range(size):
            for j in range(size):
                #If a site is occupied 
                if N_ij[i][j] == 1:
                    # --------------------------new cluster ---------------------
                    #If a site is occupied and is in corner - found new cluster
                    if i == 0 and j == 0:
                        N_ij[i][j] = k
                        M_k[k] += 1
                        k += 1

                    #If a site is occupied and the left neighbor is empty (first row) - found new cluster
                    elif i == 0 and N_ij[i][func.PI(j)] == 0:
                        N_ij[i][j] = k
                        M_k[k] += 1
                        k += 1

                    #If a site is occupied and the top neighbor is empty (first col) - found new cluster
                    elif j == 0 and N_ij[func.PI(i)][j] == 0:
                        N_ij[i][j] = k
                        M_k[k] += 1
                        k += 1

                    #If a site is occupied and the top and left neighbors are empty [row][col] - found new cluster
                    elif N_ij[func.PI(i)][j] == 0 and N_ij[i][func.PI(j)] == 0:
                        N_ij[i][j] = k
                        M_k[k] += 1
                        k += 1

                    #-------------------------- adding site to clusters --------------------
                    # top and left are different clusters
                    elif abs(M_k[N_ij[func.PI(i)][j]]) > 0 and abs(M_k[N_ij[i][func.PI(j)]]) > 0 and N_ij[i][func.PI(j)] != N_ij[func.PI(i)][j]:
                        # #choosing one of them (top)
                        N_ij[i][j] = N_ij[func.PI(i)][j]
                        k1 = N_ij[func.PI(i)][j]
                        M1 = M_k[N_ij[func.PI(i)][j]]
                        k2 = N_ij[i][func.PI(j)]
                        M2 = M_k[N_ij[i][func.PI(j)]]
                        #search for parents
                        while M1 <= 0:
                            k1 = -M1
                            M1 = M_k[k1]
                        while M2 <= 0:
                            k2 = -M2
                            M2 = M_k[k2]
                        #if left is child of top
                        if k1==k2:
                            M_k[k1]+=1
                        else:
                            #add masses
                            M_k[k1] = M_k[k1] + M_k[k2] + 1
                            #mark the second array
                            # m_k2 <- -k1 
                            M_k[k2] = -k1
                    
                    # If one of the sites (top or left) has the value of another cluster
                    # top is different
                    elif N_ij[func.PI(i)][j] > 1:
                        #choosing one of them (top)
                        N_ij[i][j] = N_ij[func.PI(i)][j]
                        k1 = N_ij[func.PI(i)][j]
                        M1 = M_k[k1]
                        #searching for parent
                        while M1 < 0:
                            k1 = -M1
                            M1 = M_k[k1]
                        #adding mass to main cluster
                        M_k[k1] += 1 

                    #left is different
                    elif N_ij[i][func.PI(j)] > 1:
                        N_ij[i][j] = N_ij[i][func.PI(j)]
                        k2 = N_ij[i][func.PI(j)]
                        M2 = M_k[k2]
                        #searching for parent
                        while M2 < 0:
                            k2 = -M2
                            M2 = M_k[k2]
                        #adding mass to main cluster
                        M_k[k2] += 1
        # M_k = np.delete(M_k, np.where(M_k < 1))
        return(N_ij,M_k)



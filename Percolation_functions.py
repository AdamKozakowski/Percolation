import Lattice
import numpy as np
import matplotlib.pyplot as plt

def MonteCarloSimulation(L,T,p0,pk,dp):
        '''Monte Carlo simulation for percolation model. It simulates T independent 
        latices of size L or range of probabilites of ocupations (poc) np.arange(p0,pk,dp).
        inputs: L: int - size of square lattice
                T: number of simulations
                p0: starting poc
                pk: ending poc
                dp: poc step
        outputs:
                output txt file  "Ave_L_{L}_T_{T}.txt" 
                with columns "p\tP_flow\tavg_s_max\n" 
                coresponding to:
                p - current poc, 
                P_flow - chance of reaching end row
                avg_s_max - avrage size of biggest cluster'''
        #creating ranges
        p_range = np.arange(p0,pk,dp)
        pflow = 0
        smax = 0
        #creating output files
        filename1 = f"Ave_L_{L}_T_{T}_pc.txt"
        with open(filename1, "w", encoding="utf8") as output_file1:
            output_file1.write("p  P_flow  avg_s_max\n")
                #for every p in range
            for p in p_range: 
                # Dictionary to accumulate cluster size counts
                cluster_sizes_aggregated = {}   
                filename2 = f"distributions/DISt_p_{round(p, 6)}_L_{L}_T_{T}_pc.txt"
                with open(filename2, "w", encoding="utf8") as output_file2:
                    #repeat T times (independent MC simulations)
                    for _ in range(T):
                        #create lattice
                        model = Lattice.Lattice(L,p)
                        #burn lattice and check if fire reached end
                        end_reached = model.burningMethod()[1]
                        if end_reached == True:
                            #increase number of reaches
                            pflow += 1
                        #Anysis of clusters sizes
                        cluster_sizes = model.HoshenKopelman()[1]
                        smax += max(cluster_sizes, default=0)
                        # Accumulate counts for each cluster size
                        for size in cluster_sizes:
                            if size in cluster_sizes_aggregated:
                                cluster_sizes_aggregated[size] += 1
                            else:
                                cluster_sizes_aggregated[size] = 1
                    #write to files
                    output_file2.write("s  n(s, p, L)\n")
                    total_cluster_sizes = sum(cluster_sizes_aggregated.values())
                    for size, count in cluster_sizes_aggregated.items():
                        output_file2.write(f"{size}  {count/total_cluster_sizes}\n")

                    #write p and its probability
                    output_file1.write(f"{p}  {pflow/T:.6f}  {smax/T:.6f}\n")

                    #restart counting
                    pflow = 0
                    smax = 0
                print(f"Analysis result saved to {filename2}")
        print(f"Analysis result saved to {filename1}")


def isEndReached(burning_forest,size, t):
    '''Function checking if connection of 1st and last row of lattice egzists
    (index of last burned site is on last row).
    output: Bool'''
    for i in range(size):
        if burning_forest[size-1,i] == t:
            return True
    return False
# next and previous indexes methods
def NI(size, i):
    '''Funtion returning index of next neighbour
    boundary condidtion - return i'''
    if i<size-1:
        return i+1
    else: 
        return i
        
def PI(i):
    '''Funtion returning index of previous neighbour
    boundary condidtion - return i'''
    if 0<i:
        return i-1
    else: 
        return i      
    
def visualize(lat,model :Lattice):
    '''Function for vizualizing modified Lattice object
    returns heatmap'''
    plt.imshow(lat, cmap='afmhot', interpolation='nearest')

    # add text annotations/labels for each element
    for i in range(lat.shape[0]):
        for j in range(lat.shape[1]):
            plt.text(j, i, str(lat[i, j]), ha='center', va='center', color='black')

    plt.colorbar()  # add colorbar for reference
    plt.title(f"L={model.size}, p={model.p}")
    plt.show()
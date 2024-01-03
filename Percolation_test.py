import Lattice
import numpy as np
import Percolation_functions as func
import matplotlib.pyplot as plt

# Initializing percolation model
model = Lattice.Lattice(10,0.67)
# func.visualize(model.lat, model)

# ---------------testing buning method:---------------
# model = Lattice.Lattice(10,0.67)
# burn, percolation, shortest_path = model.burningMethod()
# func.visualize(burn, model)
# print("Percolation reached: ", percolation)
# if percolation:
#     print("Shortest path: ", shortest_path)

# -------------testing of HK algorythm'------------------
#-----new cluster-------------
# #corner
# model.__setLat__(np.array([[1,0,0],[0,0,0],[0,0,0]]))
# #top band
# model.__setLat__(np.array([[0,1,0],[0,0,0],[0,0,0]]))
#left band
#model.__setLat__(np.array([[0,0,0],[1,0,0],[0,0,0]]))
#no neighbours
# model.__setLat__(np.array([[0,0,0],[0,1,0],[0,0,0]]))

#-----adding to clusters-----------------
# top and left are different clusters (positive masses)
# model.__setLat__(np.array([[0,1,0],[1,1,0],[0,0,0]])) #[ 0  0  3 -2]
#top and left are different clusters (negative masses)
# model.__setLat__(np.array([[1,0,1,0,1,0,1],[1,0,1,0,1,1,1],[1,1,1,0,1,0,0],[1,0,0,0,1,0,0],[1,1,1,1,1,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]])) #[ 0  0 -3 -5 -5 20  0
# top (negative mass) and left (positive mass) are different clusters 
# model.__setLat__(np.array([[1,0,1,0,1],[1,0,1,1,1],[1,1,1,0,0],[0,0,0,0,0],[0,0,0,0,0]])) #[ 0  0 -4 -4 10 
# top (positive mass) and left (negative mass) are different clusters 
#model.__setLat__(np.array([[1,0,1,0,1],[1,1,1,0,1],[1,0,0,0,1],[1,1,1,1,1],[0,0,0,0,0]])) #[ 0  0 -3 -4 14

# --------If one of the sites (top or left) has the value of another cluster ----------
# top is different (positive)
#model.__setLat__(np.array([[0,1,0],[0,1,0],[0,0,0]])) #[0 0 2 0]
# top is different (negative)
#model.__setLat__(np.array([[1,0,1],[1,1,1],[1,0,0]])) #[ 0  0 -3  6]
#left is different(positive)
#model.__setLat__(np.array([[0,0,0],[0,1,1],[0,0,0]])) #[0 0 2 0]
# left is different (negative)
#model.__setLat__(np.array([[1,0,1,0],[1,1,1,0],[1,0,0,0],[1,1,0,0]])) #[ 0  0 -3  8  0
#top and left are the same
#model.__setLat__(np.array([[0,1,1],[0,1,1],[0,0,0]])) #[0 0 4 0]
#bug - top main cluster (positive M1), left child cluster (negative mass m2=-k1)
model.__setLat__(np.array([[1,1,0,0,0,0],[1,0,1,0,1,0],[1,1,1,1,1,1],[0,0,1,1,1,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])) #[ 0  0 -5  1  1 -6 -6  0
lat2, sizes = model.HoshenKopelman()
print(sizes)
func.visualize(lat2,model)
#--------------------- b) 
# testing MC simulation 
# model = Lattice.Lattice(10,0.3)
# print(model.__getsize__)
# lat2, sizes = model.HoshenKopelman()
# print(sizes)
# func.visualize(lat2,model)

# a)
#1. check if the path connecting the first and the last row exists
#testing buning method:
# burn, percolation, shortest_path = model.burningMethod()
# print(model.burningMethod()[1])
# Lattice.visualize(burn, model)
# print("Percolation reached: ", percolation)
# if percolation:
#     print("Shortest path: ", shortest_path)


# #2. The maximum cluster size smax
# lat2, sizes = model.HoshenKopelman()
# print("Maximal cluster size in given lattice: ", max(sizes))
# Lattice.visualize(lat2,model)

# 3. distribution of clusters n(s, p, L), where s is a size of a cluster (to find clusters use the Hoshen-
# Kopelman (HK) algorithm, check e.g. pages 28-31 in ICP )
# model = Lattice.Lattice(10,0.3)
# lat2, sizes = model.HoshenKopelman()
# plt.hist(sizes, bins=model.__getsize__()**2)
# plt.xlabel("Cluster size s")
# plt.ylabel("Number of clusters")
# plt.title(f"Distribution of clusters for L = {model.__getsize__()}, p = {model.__getp__()}")
# plt.show()

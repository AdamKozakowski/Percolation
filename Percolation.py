import numpy as np
import Percolation_functions as func

# Using Monte Carlo simulations, which should consist of T trials generate the output files for:
# Probability Pflow that the path connecting the first and the last row exists as a function of p
# for L = 10, 50, 100 (use legend)

#All this data should be saved in perc−ini.txt and used to run the program (program should load this
#data as an input).
L,T,p0,pk,dp = 0, 0, 0, 0, 0
with open('perc_init.txt') as file:
    next(file) #skipping 1st line
    for line in file:
        fields = line.strip().split()
        L,T,p0,pk,dp = int(fields[0]), int(fields[1]), float(fields[2]), float(fields[3]), float(fields[4])

func.MonteCarloSimulation(L,T,p0,pk,dp)
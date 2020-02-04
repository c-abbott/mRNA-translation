from ProteinSynthesis import ProteinSynthesis
import sys
import numpy as np

def main():
    if len(sys.argv) != 2:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <parameters file>")
        quit()
    else:
        infile_parameters = sys.argv[1]

    # open input file and assinging parameters
    with open(infile_parameters, "r") as input_file:
        # read the lines of the input data file
        line = input_file.readline()
        items = line.split(", ")

        l = int(items[0])  # Ribosome length.
        L = int(items[1])  # Lattice length.
        alpha = float(items[2])   # Initiation rate.
        mcsteps = int(items[3])   # Monte Carlo steps.
    omegas = np.ones(L-1)*1.1 # Transition Rates
    omegas[-1] = 1000

    simulation = ProteinSynthesis(length = l, size = L, alpha = alpha, omegas = omegas)
    
    for i in range(mcsteps):
        R = simulation.get_R()
        print(R)
        t = simulation.get_random_time(R)
        index = simulation.get_transition(R)
        simulation.update(index)
        print(index)
        print(simulation.taus[0:15])
        print(simulation.a[0:15])
        

main()

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

        r_length = int(items[0])  # Ribosome length.
        l_length = int(items[1])  # Lattice length.
        alpha = float(items[2])   # Initiation rate.
        mcsteps = int(items[3])   # Monte Carlo steps.
    omegas = np.ones(l_length-1)*0.8 # Transition Rates

    simulation = ProteinSynthesis(length = r_length, size = l_length, alpha = alpha, t_rates = omegas)
    
    for i in range(mcsteps):
        R = simulation.get_R()
        t = simulation.get_random_time(R)
        index = simulation.get_transition(R)
        print(index)






main()

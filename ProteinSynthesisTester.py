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
        omegas = int(items[3])    # Transition rates (includes beta).

    simulation = ProteinSynthesis(length = r_length, size = l_length, alpha = alpha, t_rates = omegas)

main()
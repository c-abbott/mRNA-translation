from ProteinSynthesis import ProteinSynthesis
import sys
import numpy as np
import matplotlib.pyplot as plt


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

        l = int(items[0])        # Ribosome length.
        L = int(items[1])        # Lattice length.
        alpha = float(items[2])  # Initiation rate.
        beta = float(items[3])   # Detach rate.
        T = int(items[4])        # Upper time limit.
        n_traj = int(items[5])   # Number of trajectories.
        n_meas = int(items[6])   # Number of measurements.
    omegas = np.ones(L)*1.0      # Transition Rates.
    omegas[0] = 0                # Nothing will occupy first site.
    omegas[-1] = beta            # Detatch at final site.

    # Denisities storage.
    densities = []

    for i in range(n_traj):
        # Create new instance of the simulation for every trajectory.
        simulation = ProteinSynthesis(length=l, size=L, alpha=alpha, omegas=omegas)
        # Setting time domains.
        dt = T / n_meas
        measure_times = np.arange(0, n_meas, dt)
        ts = []
        t = 0
        while t < T:
            # Collect all possible moves.
            R = simulation.get_R()
            # Sample random time
            t += simulation.get_random_time
            ts.append(t)
            # Choose which ribosome moves.
            index = simulation.get_transition(R)
            # Update simulation - ribosome hops.
            simulation.update(index)
    

    


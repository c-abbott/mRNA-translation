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

        l = int(items[0])        # Ribosome length.
        L = int(items[1])        # Lattice length.
        alpha = float(items[2])  # Initiation rate.
        beta = float(items[3])   # Detach rate.
        mcsteps = int(items[4])  # Monte Carlo steps.
        n = int(items[5])        # MC sampling frequency.
        tol = float(items[6])    # Steady state tolerance.
    omegas = np.ones(L)*1.1      # Transition Rates.
    omegas[0] = 0                # Nothing will occupy first site.
    omegas[-1] = beta            # Detatch at final site.

    # Create instance of the simulation.
    simulation = ProteinSynthesis(length = l, size = L, alpha = alpha, omegas = omegas)
    # Setting time domain
    time_steps = np.zeros(mcsteps) # Array of all time step sizes.
    times = np.zeros(mcsteps) # Cumulative sum of times.
    # Setting densities.
    density = np.zeros(simulation.size)
    # Initially no ribosomes on mRNA strand.
    occ_num = 0
    # For determining time to reach steady state.
    outcome = False

    # Simulation begins.
    for i in range(mcsteps):
        # Collect all possible moves.
        R = simulation.get_R()
        # Increase time.
        t = simulation.get_random_time(R)
        time_steps[i] = t
        times[i] = np.sum(time_steps[:i+1])
        # Choose which ribosome moves.
        index = simulation.get_transition(R)
        # Update simulation - ribosome hops.
        simulation.update(index)

        # Determine time to reach steady state.
        if i % n == 0 and outcome == False:
            occ_num_new = simulation.get_occupation_number()
            outcome = simulation.ss_test(occ_num_new, occ_num, tol)
            occ_num = occ_num_new
            if outcome == True:
                ss_time = times[i]

main()

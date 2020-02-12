from ProteinSynthesis import ProteinSynthesis
import sys
import numpy as np
import matplotlib.pyplot as plt
import copy

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
    omegas = np.ones(L)*1.0      # Transition Rates.
    omegas[0] = 0                # Nothing will occupy first site.
    omegas[-1] = beta            # Detatch at final site.

    # Create instance of the simulation.
    simulation = ProteinSynthesis(length = l, size = L, alpha = alpha, omegas = omegas)
    # Setting time domains.
    time_steps = np.zeros(mcsteps) # Array of all time step sizes.
    times = np.zeros(mcsteps)      # Cumulative sum of times.
    # Setting densities and initial state.
    state = np.zeros(simulation.size)
    densities = np.zeros(simulation.size)
    currents = np.zeros(simulation.size)
    # Initial number of ribosomes on mRNA strand.
    occ_num = 0
    # For determining time to reach steady state.
    outcome = False
    ss_time = 0 # Steady state time.

    # Simulation begins.
    for i in range(mcsteps):
        print(i)
        # Collect all possible moves.
        R = simulation.get_R()
        # Increase time.
        time_steps[i] = simulation.get_random_time(R)
        times[i] = np.sum(time_steps[:i+1])
        # Choose which ribosome moves.
        index = simulation.get_transition(R)
        # Update simulation - ribosome hops.
        simulation.update(index)
        # Collect data after steady state reached.
        if ss_time != 0:
            densities += simulation.get_densities(state,
                                                  (times[i]-ss_time), (times[i-1]-ss_time))
            currents[index] += 1
        # Store updated state.
        state = np.array(simulation.taus)
    
        # Determine time to reach steady state.
        if i % n == 0 and outcome == False:
            occ_num_new = simulation.get_occupation_number()
            outcome = simulation.ss_test(occ_num_new, occ_num, tol) 
            occ_num = occ_num_new
            if outcome == True:
                ss_time = times[i]

    # Time spent in steady state
    t_s = times[-1]-ss_time
    print(t_s)
    print(np.mean(densities[1:]/t_s))
    print(np.mean(currents/t_s))
    # Plotting
    simulation.plot_density(np.arange(1,simulation.size,1), densities[1:]/t_s)
    simulation.plot_current(np.arange(0,simulation.size,1), currents/t_s)
main()

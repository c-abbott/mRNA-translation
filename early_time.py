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
    omegas = np.ones(L)*5.0      # Transition Rates.
    omegas[0] = 0                # Nothing will occupy first site.
    omegas[-1] = beta            # Detatch at final site.
    # Setting time domains for observables.
    dt = T / n_meas
    measure_times = np.arange(0, T, dt)
    # Data storage
    densities = []

    for i in range(n_traj):
        # Data storage.
        traj_densities = []
        # Time initiation.
        t_old = 0
        t_new = 0
        print(i)
        # Create new instance of the simulation for every trajectory.
        simulation = ProteinSynthesis(
            length=l, size=L, alpha=alpha, omegas=omegas)
        # Begin trajectory.
        while t_new <= T:
            # Collect all possible moves.
            R = simulation.get_R()
            # Sample random time.
            t_new += simulation.get_random_time(R)
            # Count time ticks.
            ticks = simulation.count_ticks((t_new - t_old), dt)
            for j in range(ticks):
                traj_densities.append(
                    simulation.get_occupation_number() / (simulation.size - 1))
            # Update t_old
            t_old = t_new
            # Choose which ribosome moves.
            index = simulation.get_transition(R)
            # Update simulation - ribosome hops.
            simulation.update(index)

        for k in range(n_meas - len(traj_densities)):
            traj_densities.append(
                simulation.get_occupation_number() / (simulation.size - 1))
        densities.append(np.array(traj_densities[:n_meas]))

    densities_array = 1 / n_traj * np.mean(densities, axis=0)
    print(densities_array)

    plt.plot(measure_times, densities_array)
    plt.show()


main()
        

        
        
    

    


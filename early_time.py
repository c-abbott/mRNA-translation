from ProteinSynthesis import ProteinSynthesis
import sys
import numpy as np
import matplotlib.pyplot as plt


def main():
    # Taking in arguments.
    if len(sys.argv) != 3:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] +
              " <parameters file>" + "<translation rates>")
        quit()
    else:
        infile_parameters = sys.argv[1]
        trans_params = sys.argv[2]

    # Open input file and assinging parameters.
    with open(infile_parameters, "r") as input_file:
        # Read the lines of the input data file.
        line = input_file.readline()
        items = line.split(", ")

        l = int(items[0])        # Ribosome length.
        alpha = float(items[1])  # Initiation rate.
        beta = float(items[2])   # Detach rate.
        T = int(items[3])        # Upper time limit.
        n_traj = int(items[4])   # Number of trajectories.
        n_meas = int(items[5])  # Number of measurements.
    
    # Initialising elongation rates.
    with open(trans_params, "r") as f:
        elong_omegas = np.array([float(line.split()[1]) for line in f])
    omegas = np.zeros(elong_omegas.size + 2)
    omegas[1:-1] = elong_omegas  # Adding elongation rates.
    omegas[0] = 0                # Nothing will occupy first site.
    omegas[-1] = beta            # Termination rate.
    L = int(omegas.size)         # Initialising length of mRNA.
    
    # Setting time domains for observables.
    dt = T / n_meas
    measure_times = np.arange(0, T, dt)
    # Data storage
    densities = []

    # Simulations begin.
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
        
        # Add on missed ticks.
        for k in range(n_meas - len(traj_densities)):
            traj_densities.append(
                simulation.get_occupation_number() / (simulation.size - 1))
            
        # Store each trajectory.
        densities.append(np.array(traj_densities[:n_meas]))

    # Compute average over all trajectories.
    densities_array = 1 / n_traj * np.mean(densities, axis=0)

    # Plotting.
    simulation.plot_density(measure_times, densities_array)

main()
        

        
        
    

    


import numpy as np
import random
import math
import matplotlib.pyplot as plt

class ProteinSynthesis(object):
    """
        A class to simulate protein synthesis. Specifically the
        process of mRNA translation.
        =======================================================
        Attributes:
        length - int, length of ribosome.
        size - int, size of mRNA strand.
        alpha - float, initiation rate of ribosome attaching to mRNA strand.
        t_rates - ndarray, array of transition rates between codon sites.
    """

    def __init__(self, length, size, alpha, omegas):
        # Initialising parameters.
        self.length = int(length)
        self.size = int(size)
        self.alpha = float(alpha)
        self.omegas = omegas  # goes from site 2 to site L (site L = beta)
        self.build_strand()
        self.build_propensity()

    def build_strand(self):
        """
            Creates empty 1D lattice acting as the 
            mRNA strand.
        """
        self.taus = np.zeros(self.size, dtype=int)

    def build_propensity(self):
        """
            Creates a propensity array in order to tell
            which state transitions are possible.
        """
        self.a = np.zeros(self.size)
        self.a[0] = self.alpha

    def get_R(self):
        """
            Calculates R - the sum of all possible 
            transitions.
        """
        R = np.sum(self.a)
        return (R)

    def get_random_time(self, R):
        """
            Generates random sample from exponential dist
            in order to perform inverse transform sampling.
        """
        return (-math.log(random.uniform(0, 1)) / R)

    def get_transition(self, R):
        """
            Finds the index to indicate which transition
            occurs i.e. which ribosome hops.
        """
        r = random.uniform(0, 1) * R
        sum = self.a[0]
        j = 0
        while sum < r:
            sum += self.a[j+1]
            j += 1
        return j
    
    def update(self, index):
        """
            Update method to realise the transition
            chosen by the Gillespie Algorithm
        """
        # Initiation.
        if index == 0:
            self.taus[1] = 1
            self.a[0] = 0
            self.a[1] = self.omegas[1]*(1 - self.taus[self.length+1])
        # Elongation at index.
        elif index >= 1 and index <= self.length:
            self.taus[index] = 0
            self.taus[index+1] = 1
            self.a[index] = 0
            self.a[index+1] = self.omegas[index+1] * \
                (1-self.taus[index+1+self.length])
        # Elongation at index and intiation.
        elif index == (self.length + 1):
            self.taus[index] = 0
            self.taus[index+1] = 1
            self.a[index] = 0
            self.a[index+1] = self.omegas[index+1] * \
                (1-self.taus[index+1+self.length])
            self.a[0] = self.alpha

        # Elongation and potential unblocking.
        elif index >= (self.length + 2) and index <= (self.size-self.length-2):
            self.taus[index] = 0
            self.taus[index+1] = 1
            self.a[index] = 0
            self.a[index+1] = self.omegas[index+1] * \
                (1-self.taus[index+1+self.length])
            self.a[index-self.length] = self.omegas[index - self.length] \
                                        * self.taus[index - self.length]

        # No ribosomes ahead.
        elif index >= (self.size-self.length-1) and index <= (self.size-2):
            self.taus[index] = 0
            self.taus[index+1] = 1
            self.a[index] = 0
            self.a[index+1] = self.omegas[index+1]
            self.a[index-self.length] = self.omegas[index - self.length] \
                                        * self.taus[index - self.length]
            
        # Detaching from lattice.
        elif index == (self.size-1):
            self.taus[index] = 0
            self.a[index] = 0
            self.a[index-self.length] = self.omegas[index - self.length] \
                                        * self.taus[index - self.length]
        
    def get_occupation_number(self):
        """
            A method to determine the total number
            of ribosomes occupying the mRNA strand.
        """
        return np.sum(self.taus)

    def ss_test(self, occ_num_new, occ_num_old, tolerance):
        """
            A method to determine the waiting time
            for the system to reach a steady state.
       """
        if (occ_num_new - occ_num_old) / occ_num_new <= tolerance:
            return True
        else:
            return False
    
    def get_densities(self, state, t1, t0):
        """
            A method to calculate the unormalised probability
            density at each lattice site.
        """
        densities = state * (t1 - t0)
        return densities
    
    
    def plot_density(self, x_data, y_data):
        """
            Density plotter for mRNA strand.
        """
        plt.title("Lattice Site Density")
        plt.xlabel("Lattice Site")
        plt.ylabel(r"Denisty [$\rho$]")
        plt.grid()
        plt.plot(x_data, y_data)
        plt.show()
    
    def plot_current(self, x_data, y_data):
        """
            Current plotter for mRNA strand.
        """
        plt.title("Lattice Site Current")
        plt.xlabel("Lattice Site")
        plt.ylabel("Current [J]")
        plt.grid()
        plt.plot(x_data, y_data)
        plt.show()

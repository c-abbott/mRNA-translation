import numpy as np
import random
import math


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

    def __init__(self, length, size, alpha, t_rates):
        # Initialising parameters.
        self.length = int(length)
        self.size = int(size)
        self.alpha = float(alpha)
        self.t_rates = t_rates  # goes from site 2 to site L (site L = beta)
        self.build_strand()
        self.build_propensity()

    def build_strand(self):
        """
            Creates empty 1D lattice acting as the 
            mRNA strand.
        """
        self.taus = np.zeros(self.size)

    def build_propensity(self):
        """
            Creates a propensity array in order to tell
            which state transitions are possible.
        """
        # Obtain initial propensity.
        #a_1 = self.alpha
        #for j in range(1, self.length + 1):
        #    a_1 *= (1 - self.taus[j])
        # Ensures dimensions are consistent.
        #taus_l = np.append(self.taus[self.length:], np.zeros(self.length - 1))
        # Array of possible transitions.
        self.a = np.zeros(self.size)
        self.a[0] = self.alpha

    def get_R(self):
        """
            Calculates R - the sum of all possible 
            transitions.
        """
        R = np.sum(self.a)
        return(R)

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
        r = random.randint(0, R)
        sum = self.a[0]
        j = 0
        while sum < r:
            sum += self.a[j]
            j += 1
        return j

    #def update_prop(self, index):
        
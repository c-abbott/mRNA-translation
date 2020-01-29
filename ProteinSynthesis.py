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
        self.length = length 
        self.size = size
        self.alpha = alpha
        self.t_rates = t_rates
        self.build_strand()

    def build_strand(self):
        """
            Creates empty 1D lattice acting as the 
            mRNA strand.
        """
        self.taus = np.zeros(self.size)

    def get_propensity(self):
        """
            Creates a propensity array in order to tell
            which state transitions are possible.
        """
        # Obtain initial propensity.
        a_1 = self.alpha
        for j in range(1, self.length+1):
            a_1 *= (1-self.taus[j])
        # Ensures dimensions are consistent.
        taus_l = np.append(self.taus[self.length:], np.zeros(self.length-1))
        # Array of possible transitions.
        a = self.t_rates*self.taus[1:]*(1 - taus_l)
        a = np.append(a_1, a)
        return a
    
    def get_R(self, a):
        """
            Calculates R - the sum of all possible 
            transitions.
        """
        R = np.sum(a)
        return(R)

    def get_random_time(self, R):
        """
            Generates random sample from exponential dist
            in order to perform inverse transform sampling.
        """
        return (-math.log(random.uniform(0,1)) / R)
    
    def get_random_trans(self, R):
        """
            Generates random number between 0 and R to
            determine which transition occurs.
        """
        return(random.uniform(0,1)*R)



       

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

        elif index >= (self.length + 2) and index <= (self.size-self.length-2):
            self.taus[index] = 0
            self.taus[index+1] = 1
            self.a[index] = 0
            self.a[index+1] = self.omegas[index+1] * \
                (1-self.taus[index+1+self.length])

        # No ribosomes ahead.
        elif index >= (self.size-self.length-1) and index <= (self.size-2):
            self.taus[index] = 0
            self.taus[index+1] = 1
            self.a[index] = 0
            self.a[index+1] = self.omegas[index+1]
            
        # Detaching from lattice.
        elif index == (self.size-1):
            self.taus[index] = 0
            self.a[index] = 0

    def unblocker(self, index):
        positions = np.where(self.taus == 1)[0]
        dist =  positions[-1] - positions[0]
        print(dist)
        if dist > self.length:
            self.a[index+1] = self.omegas[index]

            
            


        

        

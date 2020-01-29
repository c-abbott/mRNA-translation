import numpy as np

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
        self.taus = np.zeros(self.size)

    def get_propensity(self):
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
        R = np.sum(a)
        return R



       

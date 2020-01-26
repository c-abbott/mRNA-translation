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

    def sum_t_rates(self):
        # Ensures dimensions are consistent.
        taus_l = np.append(self.taus[self.length:], np.zeros(self.length))
        # Array of possible transitions.
        a = self.t_rates*self.taus*(1 - taus_l)
        # Calculates sum of allowed transition rates.
        R = np.sum(a)


       

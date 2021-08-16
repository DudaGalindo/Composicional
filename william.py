# Test file for call lnphi calculation
from equation_of_state import PengRobinson
import numpy as np

class Properties:

    def __init__(self, w, Bin, R, Tc, Pc, T):
        self.w = w
        self.Bin = Bin
        self.R = R
        self.Tc = Tc
        self.Pc = Pc
        self.T = T
        self.Nc = len(w)


kprop = Properties(w, Bin, R, Tc, Pc, T)
PR = PengRobinson(kprop)
ph = np.array([1]) # se for líquido (x), ph = 1 e se vapor (y), ph = 0
lnphi = PR.lnphi(kprop, x, P, ph)

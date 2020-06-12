import numpy as np
from relative_permeability2 import BrooksAndCorey
#import unittest

class Tests_Relative_Perm:
    def test1(self):
        krw0 = 0.15
        kro0 = 0.5
        krg0 = 0.85

        Swr = 0.3
        Sgr = 0.35
        #Sor = 0.25

        n_w = 3
        n_o = 3
        n_g = 2
        BrooksAndCorey(Sorw, Sorg, Sgr, Swc, n_w, n_o, n_g, krw0, krg0, kro0)

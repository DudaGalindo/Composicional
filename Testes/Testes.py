import numpy as np
from ..Teste_Estabilidade import StabilityTest
import unittest

class estabilidade(unittest.TestCase):

    def test_caso2(self):
        z = np.array([0.2,0.8]) #exemplo aleatório
        caso = 2

        if caso == 2:
            Nc = 2;
            P = (100 * 10E5 / 101325) * 14.7;
            R = 10.73159;
            T = 350 *(9/5);
            Tc = (9/5)*np.array([190.6, 460.4]);
            Pc = 14.7*np.array([45.4,33.4]);
            w = np.array([0.008,0.227]);
            Bin = np.array([[0,0.0236],[0.0236,0]]);

        estavel = StabilityTest.Stability(w,Bin,R,Tc,Pc,T,P,Nc,z)

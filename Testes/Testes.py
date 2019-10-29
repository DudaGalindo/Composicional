import numpy as np
from ..Teste_Estabilidade import StabilityTest
import unittest

class estabilidade(unittest.TestCase):
#Unidades estão atualmente no SI
    def test_caso2(self):
        z = np.array([0.2,0.8]) #exemplo aleatório
        caso = 2

        if caso == 2:
            Nc = 2;
            P = (100)*10E5 #* 10E5 / 101325) * 14.7; pressão de 100bar e ele converte para atm e depois psi
            R = 8.3144598 #10.73159; - unid inglesa
            T = 350 + 273#*(9/5); - T em C originalmente - faltou ele somar 32
            Tc = np.array([190.6, 460.4]) + 273; #(9/5)*
            Pc =np.array([45.4,33.4])*10E5; # 14.7*
            w = np.array([0.008,0.227]);
            Bin = np.array([[0,0.0236],[0.0236,0]]);

        StabilityTest.Stability(w,Bin,R,Tc,Pc,T,P,Nc,z)
        StabilityTest.TPD(Nc,T,P,R,Tc,Pc,Bin,w,z)

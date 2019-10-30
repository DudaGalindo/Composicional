import numpy as np
from ..Teste_Estabilidade import StabilityTest
import unittest

R = 8.3144598

class estabilidade(unittest.TestCase):
#Unidades estão atualmente no SI

    def test_caso2(self):
        z = np.array([0.2,0.8]) #exemplo aleatório
        caso = 2

        if caso == 2:
            Nc = 2;
            P = (100)*10E5 #* 10E5 / 101325) * 14.7; pressão de 100bar e ele converte para atm e depois psi
            T = 350 #*9/5 - T em K to R
            Tc = np.array([190.6, 460.4]);
            Pc =np.array([45.4,33.4])*10E5; # 14.7*
            w = np.array([0.008,0.227]);
            Bin = np.array([[0,0.0236],[0.0236,0]]);
        C7 = 'n'
        ph = 'l'
        print('caso2:')
        StabilityTest.Stability(w,Bin,R,Tc,Pc,T,P,Nc,C7,z,ph)
        #StabilityTest.TPD(Nc,C7,T,P,R,Tc,Pc,Bin,w,z)

    def teste_table16_1(self): #checando o fator de compressibilidade - Z=1.3225 (ans:1.32411) - acho aceitavel considerando as aprox
        #               Met       Prop     n-Pent       n-Dec      n-Hexadec
        z   = np.array([0.5449,  0.1394,   0.1314,     0.0869,      0.0974])
        Tc  = np.array([   343,     666,      845,       1112,        1291])*(5/9)
        Pc  = np.array([ 667.2,   615.8,    489.4,      305.7,       205.7])*6894.7573
        w   = np.array([ 0.008,   0.152,    0.251,      0.490,       0.742])
        Bin = np.array([[0,0.009,0.021,0.052,0.080],[0.009,0,0.003,0.019,0.039],\
                        [0.021,0.003,0,0.008,0.022],[0.052,0.019,0.008,0,0.004],\
                        [0.08,0.039,0.022,0.004,0]])
        Nc = len(z)
        T = 338.7
        P = 5000*6894.7573
        # Is there heptane plus in the composition: if yes, C7 = 'y', no, C7 = 'n'
        C7 = 'y'
        ph = 'l'
        print('16_1:')
        StabilityTest.Stability(w,Bin,R,Tc,Pc,T,P,Nc,C7,z,ph)


    def teste_table16_3x(self): #stable
        #               Met       n-Pent
        z  = np.array([  0.7,      0.3])
        Tc = np.array([  343,      845])*(5/9)
        Pc = np.array([667.2,    489.4])*6894.7573
        w  = np.array([0.008,    0.251])
        Bin = np.array([[0,0.021],[0.021,0]])
        Nc = len(z)

        T = (100+459.67)*5/9
        P = 2600*6894.7573
        ph = 'l'
        C7 = 'n'
        print('16_3:')
        StabilityTest.Stability(w,Bin,R,Tc,Pc,T,P,Nc,C7,z,ph)

    def teste_table16_5z(self): #unstable
        #               Met       Prop     n-Pent       n-Dec      n-Hexadec
        z   = np.array([0.8232,  0.0871,   0.0505,     0.0198,      0.0194])
        Tc  = np.array([   343,     666,      845,       1112,        1291])*(5/9)
        Pc  = np.array([ 667.2,   615.8,    489.4,      305.7,       205.7])*6894.7573
        w   = np.array([ 0.008,   0.152,    0.251,      0.490,       0.742])
        Bin = np.array([[0,0.009,0.021,0.052,0.080],[0.009,0,0.003,0.019,0.039],\
                        [0.021,0.003,0,0.008,0.022],[0.052,0.019,0.008,0,0.004],\
                        [0.08,0.039,0.022,0.004,0]])
        Nc = len(z)
        T = (100+459.67)*5/9
        P = 1500*6894.7573
        ph = 'l/g'
        C7 = 'y'
        print('16_5:')
        StabilityTest.Stability(w,Bin,R,Tc,Pc,T,P,Nc,C7,z,ph)

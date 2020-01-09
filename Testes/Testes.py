import numpy as np
from ..Checa_Estabilidade import StabilityTest
import thermo
import unittest

#R = 8.3144598
R = 10.73159
class estabilidade(unittest.TestCase):

    def test_caso1(self):
        #Units: T[Rankine];P[psia]
        T = (200 + 459.67) #Rankine*5/9=Kelvin
        P = 160 #psi*6894.7573 = Pa
        Tc = np.array([-116.32,206.24,305.96,454.1,660.38])
        Tc = (Tc + 459.67)
        Pc = np.array([667.19386,615.75821,551.09625,477.028914,367.5444663])#*6894.7573
        w  = np.array([0.008,0.152,0.193,0.27504,0.443774])
        Mw = np.array([16.043,44.097,58.124,86,134]) #molecular weight
        Bin  = np.array([[0,0,0,0.01,0.02], [0,0,0,0.01,0.01], [0,0,0,0,0.01], [0.01,0.01,0,0,0.01], [0.02,0.01,0.01,0.01,0]])
        z = np.array([0.6,0.1,0.1,0.1,0.1])
        C7 = 0

        print('\ncaso1:')

        obj = StabilityTest(w,Bin,R,Tc,Pc,T,P,C7)
        sp1,sp2 = obj.Stability(z)
        if sp1>1 or sp2>1:
            obj.molar_properties(z,Mw)

        #checar se isso abaixo iria funcionar e se é necessário.
        #Acho que só vou saber disso se der erro sem.
        '''if sp1<1 and sp2<1:
            TPD = obj.TPD(z)
            if TPD.any()<0:
                obj.molar_properties(z,Mw)'''

        print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L[1],'V: ',obj.V[1])
        print('fl: ',obj.fl,'fv: ',obj.fv)

    def test_book(self):
        z = np.array([0.8232,0.0871,0.0505,0.0198,0.0194])
        Tc = np.array([343,666,845,1112,1291])#Rankine
        Pc = np.array([667.2,615.8,489.4,305.7,205.7]) #psia
        Bin = np.array([[0,0.009,0.021,0.052,0.080],\
                        [0.009,0,0.003,0.019,0.039],\
                        [0.021,0.003,0,0.008,0.022],\
                        [0.052,0.019,0.008,0,0.004],\
                        [0.08,0.039,0.022,0.004,0]])
        w = np.array([0.008,0.152,0.251,0.490,0.742])
        Mw = np.array([16.043,44.097,72.15,142.29,226.41])
        T = 100+459.67#Fahrenheit
        P = 1500 #psia
        C7 = 0

        print('\ncasoN:')

        obj = StabilityTest(w,Bin,R,Tc,Pc,T,P,C7)
        sp1,sp2 = obj.Stability(z)
        if sp1>1 or sp2>1:
            obj.molar_properties(z,Mw)

        print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L[1],'V: ',obj.V[1])
        print('fl: ',obj.fl,'fv: ',obj.fv)

    '''def test_caso2(self):
        Bin = np.array([[0, 0.033],[0.033, 0]])
        Tc = np.array([369.8, 425.2])*9/5
        Pc = 14.7*np.array([41.9, 37.5])
        w = np.array([0.152, 0.193])
        Mw = np.array([44.097, 58.124])
        z = np.array([0.5, 0.5])
        T = 9/5*396
        P = 14.7*(3.86e6/101325)
        C7 = 0

        print('\ncaso2:')
        obj = StabilityTest(w,Bin,R,Tc,Pc,T,P,C7)
        sp1,sp2 = obj.Stability(z)
        if sp1>1 or sp2>1:
            obj.molar_properties(z,Mw)
        print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L[1],'V: ',obj.V[1])
        print('fl: ',obj.fl,'fv: ',obj.fv)

    def test_caso3(self):
        z = np.array([0.5,0.5]) #exemplo aleatório
        P = (100*1E5)/101325*14.7# pressão de 100bar e ele converte para atm e depois psi
        T = 350*9/5 #- T em K to R
        Tc = np.array([190.6, 460.4])*9/5;
        Pc =np.array([45.4,33.4])*14.7; # 14.7*
        w = np.array([0.008,0.227])
        Bin = np.array([[0,0.0236],[0.0236,0]])
        Mw = np.array([44.097,58.124])
        C7 = 0

        print('\ncaso3:')
        obj = StabilityTest(w,Bin,R,Tc,Pc,T,P,C7)
        sp1,sp2 = obj.Stability(z)
        if sp1>1 or sp2>1:
            obj.molar_properties(z,Mw)'''

    '''if sp1<1 and sp2<1:
            TPD = obj.TPD(z)
            if TPD.any()<0: #checar se isso iria funcionar
                obj.molar_properties(z,Mw)'''

    '''    print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L[1],'V: ',obj.V[1])
        print('fl: ',obj.fl,'fv: ',obj.fv)'''

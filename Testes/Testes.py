import numpy as np
from ..Checa_Estabilidade import StabilityTest
import thermo
import unittest

R = 8.3144598

class estabilidade(unittest.TestCase):
#Unidades estão atualmente no SI
    def test_caso1(self):
        T = (200 + 459.6)*5/9 #Rankine*5/9=Kelvin
        P = 160*6894.7573 #psi to Pa
        Tc = np.array([-116.32,206.24,305.96,454.1,660.38])
        Tc = (Tc + 459.67)*5/9
        Pc = np.array([667.19386,615.75821,551.09625,477.028914,367.5444663])*6894.7573
        #Vshift = [0 0 0 0 0];
        w  = np.array([0.008,0.152,0.193,0.27504,0.443774])
        Mw = np.array([16.043,44.097,58.124,86,134]) #molecular weight
        Bin  = np.array([[0,0,0,0.01,0.02], [0,0,0,0.01,0.01], [0,0,0,0,0.01], [0.01,0.01,0,0,0.01], [0.02,0.01,0.01,0.01,0]])
        z = np.array([0.6,0.1,0.1,0.1,0.1])
        C7 = 0
        print('caso1:')

        obj = StabilityTest(w,Bin,R,Tc,Pc,T,P,C7)
        sp1,sp2 = obj.Stability(z)
        if sp1>1 or sp2>1:
            obj.molar_properties(z,Mw)

        '''if sp1<1 and sp2<1:
            TPD = obj.TPD(z)
            if TPD.any()<0: #checar se isso iria funcionar
                obj.molar_properties(z,Mw)'''

        print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L[1],'V: ',obj.V[1])
        print('fl: ',obj.fl,'fv: ',obj.fv)

    '''def test_caso2(self):
        z = np.array([0.5,0.5]) #exemplo aleatório
        P = (100*1E5)#/101325)*14.7# pressão de 100bar e ele converte para atm e depois psi
        T = 350#*9/5 #- T em K to R
        Tc = np.array([190.6, 460.4])#*9/5;
        Pc =np.array([45.4,33.4])*1E5#*14.7; # 14.7*
        w = np.array([0.008,0.227])
        Bin = np.array([[0,0.0236],[0.0236,0]])
        Mw = np.array([44.097,58.124])
        C7 = 0
        print('caso2:')
        obj = StabilityTest(w,Bin,R,Tc,Pc,T,P,C7)
        sp1,sp2 = obj.Stability(z)
        if sp1>1 or sp2>1:
            obj.molar_properties(z,Mw)'''

    '''if sp1<1 and sp2<1:
            TPD = obj.TPD(z)
            if TPD.any()<0: #checar se isso iria funcionar
                obj.molar_properties(z,Mw)'''

    '''print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L[1],'V: ',obj.V[1])
        print('fl: ',obj.fl,'fv: ',obj.fv)'''

import numpy as np
from ..Checa_Estabilidade import StabilityTest
import thermo
import unittest

R = 8.3144598

class estabilidade(unittest.TestCase):
#Unidades estão atualmente no SI
# Esses testes devem conter apenas uma fase. Se o sistema for composto por mais de uma fase, deve-se tomar x e y separadamente
    def test_caso1(self):
        T = (200 + 459.6)*5/9 #Rankine*5/9=Kelvin
        P = 160*6894.7573 #psi to Pa
        Nc = 5;
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

        obj = StabilityTest(w,Bin,R,Tc,Pc,T,P,Nc,C7)
        sp1,sp2 = obj.Stability(z)
        if sp1>1 or sp2>1:
            obj.molar_properties(z,Mw)

    def test_caso2(self):
        z = np.array([0.5,0.5]) #exemplo aleatório
        Nc = 2;
        P = (100*1E5)#/101325)*14.7# pressão de 100bar e ele converte para atm e depois psi
        T = 350#*9/5 #- T em K to R
        Tc = np.array([190.6, 460.4])#*9/5;
        Pc =np.array([45.4,33.4])*1E5#*14.7; # 14.7*
        w = np.array([0.008,0.227])
        Bin = np.array([[0,0.0236],[0.0236,0]])
        Mw = np.array([44.097,58.124])
        C7 = 0
        ph = 'l'
        print('caso2:')
        obj = StabilityTest(w,Bin,R,Tc,Pc,T,P,Nc,C7)
        sp1,sp2 = obj.Stability(z)
        if sp1>1 or sp2>1:
            obj.molar_properties(z,Mw)
        if sp1<1 and sp2<1:
            TPD = obj.TPD(z)
            if TPD.any()<0: #checar se isso iria funcionar
                obj.molar_properties(z,Mw)
        #print('x: ','y: ',obj.x,obj.y)
        #print('K: ',obj.K)
        #print('L: ','V: ',obj.L,obj.V)
        #print('fl: ','fv: ',obj.fl,obj.fv)

    '''def teste_table16_1(self): #checando o fator de compressibilidade - Z=1.3225 (ans:1.32411) - acho aceitavel considerando as aprox
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
        C7 = 1
        print('16_1:')
        obj = StabilityTest(w,Bin,R,Tc,Pc,T,P,Nc,C7,z)
        StabilityTest.run(obj)
        #StabilityTest.Stability(w,Bin,R,Tc,Pc,T,P,Nc,C7,z)


    def teste_table16_3x(self):
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
        C7 = 0
        print('16_3:')
        obj = StabilityTest(w,Bin,R,Tc,Pc,T,P,Nc,C7,z)
        StabilityTest.run(obj)

    def teste_table16_7x(self): #understand the meaning of having a mixture as the entry parameter
        #               Met       Prop     n-Pent       n-Dec      n-Hexadec
        x   = np.array([0.35630,0.14651,  0.22041,     0.13818,    0.13860])
        y   = np.array([0.89911,0.07744,  0.02288,     0.00055,    0.00002])
        nl = 0.1398;nv = 0.8602; z = x*nv+y*nl
        Tc  = np.array([   343,     666,      845,       1112,        1291])*(5/9)
        Pc  = np.array([ 667.2,   615.8,    489.4,      305.7,       205.7])*6894.7573
        w   = np.array([ 0.008,   0.152,    0.251,      0.490,       0.742])
        Bin = np.array([[0,0.009,0.021,0.052,0.080],[0.009,0,0.003,0.019,0.039],\
                        [0.021,0.003,0,0.008,0.022],[0.052,0.019,0.008,0,0.004],\
                        [0.08,0.039,0.022,0.004,0]])
        Nc = len(z)
        T = (100+459.67)*5/9
        P = 1500*6894.7573
        C7 = 1
        #print('16_7x:')
        #eos = thermo.eos_mix.SRKMIX(T=T,P=P,Tcs=Tc,Pcs=Pc,omegas=w,zs=z,kijs=Bin)
        #eosl = thermo.eos_mix.PRMIX(Tcs=Tc,Pcs=Pc,omegas=w,zs=x,kijs=Bin,T=T,P=P)
        #print('x',eosl.phase)
        #StabilityTest.Stability(w,Bin,R,Tc,Pc,T,P,Nc,C7,y)
        #eosg = thermo.eos_mix.PRMIX(Tcs=Tc,Pcs=Pc,omegas=w,zs=y,kijs=Bin,T=T,P=P)
        #print('y',eosg.phase)'''

    '''I've notice that if the mixture is in the two phase region, the EOS
    model does not identify that, in fact it only works for one phase at a time.
     So, due to that, I conclude that for this to work, the entry parameters have to
    be for each phase, like this example shows, and the stability test is made
    separetedly.'''
    '''Furthermore I'll implement the code to obtain the two phases composition
    based on having the general one (z)'''

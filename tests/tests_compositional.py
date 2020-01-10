import numpy as np
from stability_check import StabilityCheck
import thermo
import unittest


class Li_et_al_table4:
    def __init__(self):
        self.R = 10.73159
        self.Tc = np.array([-116.59, 305.69, 453.65, 652.01])
        self.Tc = (self.Tc + 459.67)
        self.Pc = np.array([667.2,551.1,430.59,305.68])
        self.w = np.array([0.008,0.193,0.296,0.49])
        self.Mw = np.array([16.043,58.124,86.16,142.29])
        self.P = 600
        self.T = 160 + 459.67
        self.Bin = np.zeros([len(self.w),len(self.w)])
        self.C7 = np.array([0,0,0,1])


class testes_Li_et_al_table4(unittest.TestCase):

    #Units: English system
    # methane n-butane n-hexane n-decane

    def test5(self):
        prop = Li_et_al_table4()
        z = np.array([-0.58,0.38,0.6,0.6])

        '''Run:'''
        obj = StabilityCheck(prop.w,prop.Bin,prop.R,prop.Tc,prop.Pc,prop.T,prop.P,prop.C7)
        obj.run(z,prop.Mw)
        
        '''Verify Results'''
        x = [0.1676,0.2239,0.3088,0.2997]
        y = [0.9115,0.0686,0.0190,0.0009]
        K = [5.4393,0.3062,0.0616,0.0031]

        for i in range(obj.Nc):
            self.assertAlmostEqual(obj.x[i],x[i],4,'ValueError:Failed')
            self.assertAlmostEqual(obj.y[i],y[i],4,'ValueError:Failed')
            if obj.x[i] != 0 and obj.y[i] != 0:
                self.assertAlmostEqual(obj.K[i],K[i],4,'ValueError:Failed')

    def test4(self):
        prop = Li_et_al_table4()
        z = np.array([0,0,0,1])

        '''Run:'''
        obj = StabilityCheck(prop.w,prop.Bin,prop.R,prop.Tc,prop.Pc,prop.T,prop.P,prop.C7)
        obj.run(z,prop.Mw)

        '''Verify Results'''
        x = [0.1724,0.0000,0.0000,0.8276]
        y = [0.9980,0.0000,0.0000,0.0020]
        K = [5.7895,0.2867,0.0538,0.0024]

        for i in range(obj.Nc):
            self.assertAlmostEqual(obj.x[i],x[i],4,'ValueError:Failed')
            self.assertAlmostEqual(obj.y[i],y[i],4,'ValueError:Failed')
            if obj.x[i] != 0 and obj.y[i] != 0:
                self.assertAlmostEqual(obj.K[i],K[i],4,'ValueError:Failed')

    def test3(self):
        prop = Li_et_al_table4()
        z = np.array([0,0,1,0])

        '''Run:'''
        obj = StabilityCheck(prop.w,prop.Bin,prop.R,prop.Tc,prop.Pc,prop.T,prop.P,prop.C7)
        obj.run(z,prop.Mw)

        '''Verify Results'''
        x = [0.1789,0.0000,0.8211,0.0000]
        y = [0.9507,0.0000,0.0493,0.0000]
        K = [5.3130,0.2999,0.0600,0.0030]

        for i in range(obj.Nc):
            self.assertAlmostEqual(obj.x[i],x[i],4,'ValueError:Failed')
            self.assertAlmostEqual(obj.y[i],y[i],4,'ValueError:Failed')
            if obj.x[i] != 0 and obj.y[i] != 0:
                self.assertAlmostEqual(obj.K[i],K[i],4,'ValueError:Failed')

    def test2(self):
        prop = Li_et_al_table4()
        z = np.array([0,1,0,0])

        '''Run:'''
        obj = StabilityCheck(prop.w,prop.Bin,prop.R,prop.Tc,prop.Pc,prop.T,prop.P,prop.C7)
        obj.run(z,prop.Mw)

        '''Verify Results'''
        x = [0.1639,0.8361,0.0000,0.0000]
        y = [0.6960,0.3040,0.0000,0.0000]
        K = [4.2473,0.3636,0.0905,0.0066]

        for i in range(obj.Nc):
            self.assertAlmostEqual(obj.x[i],x[i],4,'ValueError:Failed')
            self.assertAlmostEqual(obj.y[i],y[i],4,'ValueError:Failed')
            if obj.x[i] != 0 and obj.y[i] != 0:
                self.assertAlmostEqual(obj.K[i],K[i],4,'ValueError:Failed')

    def test1(self):
        prop = Li_et_al_table4()
        z = np.array([1,0,0,0])

        '''Run:'''
        obj = StabilityCheck(prop.w,prop.Bin,prop.R,prop.Tc,prop.Pc,prop.T,prop.P,prop.C7)
        obj.run(z,prop.Mw)

        '''Verify Results'''
        x = [0.1724,0.0000,0.0000,0.8276]
        y = [0.9980,0.0000,0.0000,0.0020]
        K = [5.7895,0.2867,0.0538,0.0024]

        for i in range(obj.Nc):
            self.assertAlmostEqual(obj.x[i],x[i],4,'ValueError:Failed')
            self.assertAlmostEqual(obj.y[i],y[i],4,'ValueError:Failed')
            if obj.x[i] != 0 and obj.y[i] != 0:
                self.assertAlmostEqual(obj.K[i],K[i],4,'ValueError:Failed')

class testes_casos_Schmall(unittest.TestCase):
    def test_caso1(self):
        #Units: T[Rankine];P[psia]
        R = 10.73159
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

        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P,C7)
        obj.run(z,Mw)

        print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L[1],'V: ',obj.V[1])
        print('fl: ',obj.fl,'fv: ',obj.fv)

    def test_caso2(self):
        R = 8.3144598
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
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P,C7)
        obj.run(z,Mw)

        print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L[1],'V: ',obj.V[1])
        print('fl: ',obj.fl,'fv: ',obj.fv)

    def test_caso3(self):
        R = 8.3144598
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
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P,C7)
        obj.run(z,Mw)

        '''if sp1<1 and sp2<1:
            TPD = obj.TPD(z)
            if TPD.any()<0: #checar se isso iria funcionar
                obj.molar_properties(z,Mw)'''

        print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L[1],'V: ',obj.V[1])
        print('fl: ',obj.fl,'fv: ',obj.fv)


    # def test_bookDandekar(self):
    #     R = 10.73159
    #     z = np.array([0.8232,0.0871,0.0505,0.0198,0.0194])
    #     Tc = np.array([343,666,845,1112,1291])#Rankine
    #     Pc = np.array([667.2,615.8,489.4,305.7,205.7]) #psia
    #     Bin = np.array([[0,0.009,0.021,0.052,0.080],\
    #                     [0.009,0,0.003,0.019,0.039],\
    #                     [0.021,0.003,0,0.008,0.022],\
    #                     [0.052,0.019,0.008,0,0.004],\
    #                     [0.08,0.039,0.022,0.004,0]])
    #     w = np.array([0.008,0.152,0.251,0.490,0.742])
    #     Mw = np.array([16.043,44.097,72.15,142.29,226.41])
    #     T = 100+459.67#Fahrenheit
    #     P = 1500 #psia
    #     C7 = np.array([0,0,0,1,1])
    #
    #     print('\ncasoN:')
    #
    #     obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P,C7)
    #     sp1,sp2 = obj.Stability(z)
    #     if sp1>1 or sp2>1:
    #         obj.molar_properties(z,Mw)
    #
    #     print('x: ',obj.x,'y: ',obj.y)
    #     print('K: ',obj.K)
    #     print('L: ',obj.L[1],'V: ',obj.V[1])
    #     print('fl: ',obj.fl,'fv: ',obj.fv)
    #
    # def test_Abbot(self):
    #     #units: SI
    #     #elements: methane ethane propane
    #     R = 8.3144598
    #     z = np.array([0.1,0.2,0.7])
    #     Tc = np.array([190.55,305.33,369.522])
    #     Pc = np.array([45.95, 48.72, 42.49])
    #     P = 13.8
    #     T = 283.15
    #     Mw = np.array([16.04,30.07,44.1])
    #     C7 = 0
    #     w = np.array([0.011,0.099,0.153])
    #     Bin = np.zeros([3,3])
    #     print('\nExemplo 10.5:')
    #
    #     obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P,C7)
    #     sp1,sp2 = obj.Stability(z)
    #     if sp1>1 or sp2>1:
    #         obj.molar_properties(z,Mw)
    #
    #     print('x: ',obj.x,'y: ',obj.y)
    #     print('K: ',obj.K)
    #     print('L: ',obj.L[1],'V: ',obj.V[1])
    #     print('fl: ',obj.fl,'fv: ',obj.fv)

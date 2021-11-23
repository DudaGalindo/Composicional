import numpy as np
from stability_check import StabilityCheck
import unittest


class Li_et_al_table4:
    @unittest.skip("ok")
    def __init__(self):
        self.R = 10.73159
        self.Tc = np.array([-116.59, 305.69, 453.65, 652.01])
        self.Tc = (self.Tc + 459.67)
        self.Pc = np.array([667.2,551.1,430.59,305.68])
        self.w = np.array([0.008,0.193,0.296,0.49])
        self.Mw = np.array([16.043,58.124,86.16,142.29])
        self.P = 600 * np.ones(5)
        self.T = 160 + 459.67
        self.Bin = np.zeros([len(self.w),len(self.w)])


class testes_Li_et_al_table4(unittest.TestCase):

    #Units: English system
    # methane n-butane n-hexane n-decane
    @unittest.skip("ok")
    def test_all(self):
        prop = Li_et_al_table4()
        z = np.array([[-0.58,0.38,0.6,0.6], [0,0,0,1], [0,1,0,0], [0,0,1,0], [1,0,0,0]])
        obj = StabilityCheck(prop.w,prop.Bin,prop.R,prop.Tc,prop.Pc,prop.T,prop.P)
        obj.run(z.T,prop.Mw)
        print('x: ', 'y: ', obj.x, obj.y)
        print('L: ','V: ', obj.L, obj.V)
        import pdb; pdb.set_trace()

    @unittest.skip("ok")
    def test5(self):
        prop = Li_et_al_table4()
        z = np.array([-0.58,0.38,0.6,0.6])[:,np.newaxis]*np.ones([len(prop.w),len(prop.P)])

        '''Run:'''
        obj = StabilityCheck(prop.w,prop.Bin,prop.R,prop.Tc,prop.Pc,prop.T,prop.P)
        obj.run(z.T,prop.Mw)

        '''Verify Results'''
        x = [0.1676,0.2239,0.3088,0.2997]
        y = [0.9115,0.0686,0.0190,0.0009]
        K = [5.4393,0.3062,0.0616,0.0031]

        for i in range(obj.Nc):
            self.assertAlmostEqual(obj.x[0,i],x[i],4,'ValueError:Failed')
            self.assertAlmostEqual(obj.y[0,i],y[i],4,'ValueError:Failed')
            if obj.x[0,i] != 0 and obj.y[0,i] != 0:
                self.assertAlmostEqual(obj.K[0,i],K[i],4,'ValueError:Failed')

    @unittest.skip("ok")
    def test4(self):
        prop = Li_et_al_table4()
        z = np.array([0,0,0,1])[:,np.newaxis]*np.ones([len(prop.w),len(prop.P)])

        '''Run:'''
        obj = StabilityCheck(prop.w,prop.Bin,prop.R,prop.Tc,prop.Pc,prop.T,prop.P)
        obj.run(z.T,prop.Mw)
        import pdb; pdb.set_trace()
        '''Verify Results'''
        x = [0.1724,0.0000,0.0000,0.8276]
        y = [0.9980,0.0000,0.0000,0.0020]
        K = [5.7895,0.2867,0.0538,0.0024]

        for i in range(obj.Nc):
            self.assertAlmostEqual(obj.x[0,i],x[i],4,'ValueError:Failed')
            self.assertAlmostEqual(obj.y[0,i],y[i],4,'ValueError:Failed')
            if obj.x[0,i] != 0 and obj.y[0,i] != 0:
                self.assertAlmostEqual(obj.K[0,i],K[i],4,'ValueError:Failed')
    @unittest.skip("ok")
    def test3(self):
        prop = Li_et_al_table4()
        z = np.array([0,0,1,0])[:,np.newaxis]*np.ones([len(prop.w),len(prop.P)])

        '''Run:'''
        obj = StabilityCheck(prop.w,prop.Bin,prop.R,prop.Tc,prop.Pc,prop.T,prop.P)
        obj.run(z.T,prop.Mw)

        '''Verify Results'''
        x = [0.1789,0.0000,0.8211,0.0000]
        y = [0.9507,0.0000,0.0493,0.0000]
        K = [5.3130,0.2999,0.0600,0.0030]

        for i in range(obj.Nc):
            self.assertAlmostEqual(obj.x[0,i],x[i],4,'ValueError:Failed')
            self.assertAlmostEqual(obj.y[0,i],y[i],4,'ValueError:Failed')
            if obj.x[0,i] != 0 and obj.y[0,i] != 0:
                self.assertAlmostEqual(obj.K[0,i],K[i],4,'ValueError:Failed')
    @unittest.skip("ok")
    def test2(self):
        prop = Li_et_al_table4()
        z = np.array([0,1,0,0])[:,np.newaxis]*np.ones([len(prop.w),len(prop.P)])

        '''Run:'''
        obj = StabilityCheck(prop.w,prop.Bin,prop.R,prop.Tc,prop.Pc,prop.T,prop.P)
        obj.run(z.T,prop.Mw)

        '''Verify Results'''
        x = [0.1639,0.8361,0.0000,0.0000]
        y = [0.6960,0.3040,0.0000,0.0000]
        K = [4.2473,0.3636,0.0905,0.0066]

        for i in range(obj.Nc):
            self.assertAlmostEqual(obj.x[0,i],x[i],4,'ValueError:Failed')
            self.assertAlmostEqual(obj.y[0,i],y[i],4,'ValueError:Failed')
            if obj.x[0,i] != 0 and obj.y[0,i] != 0:
                self.assertAlmostEqual(obj.K[0,i],K[i],4,'ValueError:Failed')

    @unittest.skip("ok")
    def test1(self):
        prop = Li_et_al_table4()
        z = np.array([1,0,0,0])[:,np.newaxis]*np.ones([len(prop.w),len(prop.P)])

        '''Run:'''
        obj = StabilityCheck(prop.w,prop.Bin,prop.R,prop.Tc,prop.Pc,prop.T,prop.P)
        obj.run(z.T,prop.Mw)

        '''Verify Results'''
        x = [0.1724,0.0000,0.0000,0.8276]
        y = [0.9980,0.0000,0.0000,0.0020]
        K = [5.7895,0.2867,0.0538,0.0024]

        for i in range(obj.Nc):
            self.assertAlmostEqual(obj.x[0,i],x[i],4,'ValueError:Failed')
            self.assertAlmostEqual(obj.y[0,i],y[i],4,'ValueError:Failed')
            if obj.x[0,i] != 0 and obj.y[0,i] != 0:
                self.assertAlmostEqual(obj.K[0,i],K[i],4,'ValueError:Failed')

class testes_casos_Schmall(unittest.TestCase):

    @unittest.skip("ok")
    def teste_caso5(self):
        R = 8.3144598
        z = np.array([1.])[:,np.newaxis]#*np.ones([10]) #exemplo aleatório
        P = np.array([ 13.78951458E6])
        #P = np.array([14619094.696946101])*np.ones(10)
        #T = np.array([15 + 273.15])
        T = np.array([366.4833])
        Tc = np.array([619.28])
        Pc =np.array([2109795.64])
        w = np.array([0.489])
        Bin = np.array([[0]])
        Mw = np.array([142.28e-3])
        C7 = np.array([0])#*np.ones([1,10])
        print('\ncaso5:')
        #import pdb; pdb.set_trace()#
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pv=0)
        obj.run(z,Mw)

        print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L,'V: ',obj.V)
        #print('fl: ',obj.fl,'fv: ',obj.fv)

    @unittest.skip("ok")
    def test1_Firoozabadi(self):
        #methane and propane
        R = 8.3144598
        z = np.array([0.3, 0.7])[:,np.newaxis] * np.ones((2,100))
        Tc = np.array([190.55, 369.522])
        Pc = np.array([4599200, 42.6e5])
        P = np.linspace(200,1200,100) * 6894.757
        T = np.array([54.4 + 273.15])
        Mw = np.array([16.04, 44.1])*1e-3
        w = np.array([0.022, 0.153])
        Bin = np.array([[0, 0.], [0., 0.]])
        #Pv = np.array([4.e6,18.97e5])
        Pv = np.array([41368542,18.97e5])
        print('\nExemplo Firoozabadi:')
        #import pdb; pdb.set_trace()
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pv)
        obj.run(z,Mw)
        import pdb; pdb.set_trace()

    @unittest.skip("?")
    def test_Firoozabadi(self):
        #methane, n-butane and n-decane
        R = 8.3144598
        z = np.array([0.7, 0.235, 0.065])[:,np.newaxis] * np.ones((3,100))
        Tc = np.array([190.55, 425.5, 619.28])
        Pc = np.array([4599200, 3.8e6, 2109795.64])
        P = np.linspace(1000,6000,100) * 6894.757
        T = np.array([71.1 + 273.15])
        Mw = np.array([16.04, 58.12, 142.28])*1e-3
        w = np.array([0.022, 0.199, 0.489])
        Bin = np.zeros((3,3))
        #Pv = np.array([4.e6,18.97e5])
        Pv = np.array([0,41368542,18.97e5])
        print('\nExemplo Firoozabadi:')
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pv)
        obj.run(z,Mw)

    @unittest.skip('ok')
    def teste_res_fluid_brb(self):
        R = 8.3144598
        #z = np.array([0.63840871, 0.09640664, 0.265184659])[:,np.newaxis]#*np.ones([3,100])
        z = np.array([0.01, 0.19, 0.8])[:,np.newaxis]*np.ones([3,100])
        Tc = np.array([304.21, 190.6, 734.68])
        Pc = np.array([7.39e6, 4.6e6, 1.74e6])
        P = np.linspace(20650000, 24700000, 100)
        #P = np.array([21279410.17965427])
        T = np.array([299.82])
        Mw = np.array([44.01e-3,16.04e-3, 222e-3])
        w = np.array([0.225, 0.0225, 0.684])
        Bin = np.array([[0, 0.12, 0.12], [0.12, 0, 0.], [0.12, 0., 0.]])
        Pv = np.array([0, 341.4e3, 0])
        Pb_guess = 9.65e6
        print('\nExemplo BRB:')
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pv)
        obj.run(z,Mw)
        import pdb; pdb.set_trace()

    @unittest.skip('ok')
    def teste_inj_fluid_brb(self):
        R = 8.3144598
        z = np.array([0.95, 0.05])[:,np.newaxis]
        Tc = np.array([304.21, 190.6])
        Pc = np.array([7.39e6, 4.6e6])
        P = np.array([101325])
        #P = np.array([20.65e6])
        Pv = np.array([0,341.4e3])
        #T = np.array([273.15])
        T = np.array([288.7])
        Mw = np.array([44.01e-3,16.04e-3])
        w = np.array([0.225,0.0225])
        Bin = np.array([[0.0,0.12],[0.12,0.0]])
        print('\nExemplo BRB:')

        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pv)
        obj.run(z,Mw)
        import pdb; pdb.set_trace()
        print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L,'V: ',obj.V)
        print('fl: ',obj.fl,'fv: ',obj.fv)
        print('ksi_L', obj.ksi_L, 'ksi_V', obj.ksi_V)

    @unittest.skip('ok')
    def test_inj_agua(self):
        R = 8.3144598
        z = np.array([0.1, 0.4, 0.4, 0.1])[:,np.newaxis]*np.ones((4,100))
        Tc = np.array([190.6, 305.4, 425.2, 543.2])
        Pc = np.array([313021.9678, 332327.2874, 258553.3874999999, 213530.62428999998])
        Mw =np.array([16.04e-3, 30.07e-3, 58.12e-3, 96e-3])
        w = np.array([0.008, 0.098, 0.193, 0.308])
        Bin = np.array([[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]])
        Pv = np.array([0.0, 0.0, 0.0, 0.0])
        P = np.linspace(10e5, 15e6, 100)
        #P = np.array([101325])
        T = np.array([288.7])
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pv)
        obj.run(z,Mw)

    @unittest.skip('ok')
    def test_water_inj_schmall3d(self):
        R = 8.3144598
        z = np.array([0.49947493593730447, 0.030015800353206768, 0.07007056389862372, 0.20021841113835698, 0.15016518896172965, 0.050055099710778636])[:,np.newaxis]
        #z = np.array([0.5, 0.03, 0.07, 0.2, 0.15, 0.05])[:,np.newaxis] #* np.ones((6,100))
        Tc = np.array([190.6, 369.8, 507.4, 617.6, 708, 768])
        Pc = np.array([4600155, 4245517.5, 2968822.5, 2107560.0, 1.47e6, 1.17e6])
        vc = np.array([0.000099, 0.000203, 0.00037, 0.000603, 0.000895, 0.00169])
        Mw = np.array([16.043e-3, 44.097e-3, 86.178e-3, 142.286e-3, 212.41e-3, 282.5e-3])
        w = np.array([0.008, 0.152, 0.299, 0.489, 0.685, 0.912])
        Bin = np.array([[0.,0.,.0,0.,0.,0.], [0.,0.,.0,0.,0.,0.], [0.,0.,.0,0.,0.,0.], [0.,0.,.0,0.,0.,0.],
                        [0.,0.,.0,0.,0.,0.], [0.,0.,.0,0.,0.,0.]])
        Pv = np.array([8e6, 0., 0., 0., 0., 0.])
        #P = np.linspace(8.46e6, 62e6, 100)
        P = np.array([206832942.5900688])
        T = np.array([344.25])
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pv)
        obj.run(z,Mw)
        import pdb; pdb.set_trace()

    @unittest.skip("ok")
    def testeeeee(self):
        R = 8.3144598
        z = np.array([0.684, 0.316])[:,np.newaxis]
        Tc = np.array([190.56, 425.16])
        Pc = np.array([4605697.676, 3799011.107])
        P = np.array([101325])
        Pv = np.array([0.0,0.0])
        T = np.array([288.7])
        Mw = np.array([16.042e-3, 58.12e-3])
        w = np.array([0.008, 0.193])
        Bin = np.array([[0., 0.03], [0.03, 0.]])
        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess)
        obj.run(z,Mw)
        import pdb; pdb.set_trace()

    @unittest.skip("ok")
    def test_caso1(self):
        #Units: T[Rankine];P[psia]
        R = 10.73159
        #R = 8.3144598
        T = (200 + 459.67)#*5/9 #Rankine*5/9=Kelvin
        P = np.array([160, 300]) #*6894.7573  #psi*6894.7573 = Pa
        Tc = np.array([-116.32,206.24,305.96,454.1,660.38])
        Tc = (Tc + 459.67)#*5/9
        Pc = np.array([667.19386,615.75821,551.09625,477.028914,367.5444663])#*6894.7573
        w  = np.array([0.008,0.152,0.193,0.27504,0.443774])
        Mw = np.array([16.043,44.097,58.124,86,134])*1e-3 #molecular weight
        Bin  = np.array([[0,0,0,0.01,0.02], [0,0,0,0.01,0.01], [0,0,0,0,0.01], [0.01,0.01,0,0,0.01], [0.02,0.01,0.01,0.01,0]])
        z = np.array([[0.6,0.1,0.1,0.1,0.1],[0.8232,0.0871,0.0505,0.0198,0.0194]])


        print('\ncaso1:')

        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P)
        obj.run(z,Mw)

        print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L,'V: ',obj.V)
        print('fl: ',obj.fl,'fv: ',obj.fv)
        import pdb; pdb.set_trace()

    @unittest.skip("ok")
    def test_caso2(self):
        R = 10.73159
        Bin = np.array([[0, 0.033],[0.033, 0]])
        Tc = np.array([369.8, 425.2])*9/5
        Pc = 14.7*np.array([41.9, 37.5])
        w = np.array([0.152, 0.193])
        Mw = np.array([44.097, 58.124])
        z = np.array([0.5, 0.5]) [:,np.newaxis]*np.ones([2,10])
        T = 9/5*396
        P = 14.7*(3.86e6/101325) * np.ones(10)


        print('\ncaso2:')
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P)
        obj.run(z.T,Mw)
        import pdb; pdb.set_trace()
        print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L,'V: ',obj.V)
        print('fl: ',obj.fl,'fv: ',obj.fv)

    @unittest.skip("ok")
    def test_caso3(self):
        R = 10.73159
        z = np.array([0.5,0.5])[:,np.newaxis]*np.ones([2,10]) #exemplo aleatório
        P = (100*1E5)/101325*14.7 * np.ones(10)# pressão de 100bar e ele converte para atm e depois psi
        T = 350*9/5 #- T em K to R
        Tc = np.array([190.6, 460.4])*9/5;
        Pc =np.array([45.4,33.4])*14.7; # 14.7*class
        w = np.array([0.008,0.227])
        Bin = np.array([[0,0.0236],[0.0236,0]])
        Mw = np.array([44.097,58.124])


        print('\ncaso3:')

        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P)
        obj.run(z.T,Mw)

        print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L,'V: ',obj.V)
        print('fl: ',obj.fl,'fv: ',obj.fv)

    @unittest.skip("not now")
    def test_bookDandekar(self):
         R = 10.73159
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
         T = np.array([100+459.67])#Fahrenheit
         P = np.array([1500]) #psia
         C7 = np.array([0,0,0,1,1])

         print('\ncasoN:')
         #import pdb; pdb.set_trace()

         obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P,C7)
         sp1,sp2 = obj.Stability(z)
         if sp1>1 or sp2>1:
             obj.molar_properties(z)

         print('x: ',obj.x,'y: ',obj.y)
         print('K: ',obj.K)
         print('L: ',obj.L[1],'V: ',obj.V[1])
         print('fl: ',obj.fl,'fv: ',obj.fv)
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
    #
    #     w = np.array([0.011,0.099,0.153])
    #     Bin = np.zeros([3,3])
    #     print('\nExemplo 10.5:')
    #
    #     obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P,C7)
    #     sp1,sp2 = obj.Stability(z)
    #     if sp1>1 or sp2>1:
    #         obj.molar_properties(z)
    #
    #     print('x: ',obj.x,'y: ',obj.y)
    #     print('K: ',obj.K)
    #     print('L: ',obj.L[1],'V: ',obj.V[1])
    #     print('fl: ',obj.fl,'fv: ',obj.fv)

class Testes_IGOR(unittest.TestCase):

    @unittest.skip("ok")
    def test_Y8(self):
        R = 8.3144598
        # C1, C2, C3, n-C5, n-C7, n-C10
        z = np.array([0.8097, 0.0566, 0.0306, 0.0457, 0.0330, 0.0244])[:,np.newaxis]
        Tc = np.array([190.6, 305.4, 369.8, 469.6, 540.3, 617.9])
        Pc = np.array([45.4, 48.2, 41.9, 33.3, 27.4, 21.0])*101325
        P = np.array([194.8])*101325
        Pv = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        T = np.array([375.3])
        Mw = np.array([16.043, 30.070, 44.097, 72.15, 100.205, 142.29])*1e-3
        w = np.array([0.008, 0.098, 0.152, 0.251, 0.305, 0.484])

        Bin = np.zeros([len(z), len(z)])

        CP1 = np.zeros_like(Tc)
        CP2 = np.zeros_like(Tc)
        CP3 = np.zeros_like(Tc)
        CP4 = np.zeros_like(Tc)

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)

        obj.run(z,Mw)

        ''' Using Thermo '''
        from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVL
        from thermo.interaction_parameters import IPDB
        import time
        constants, properties = ChemicalConstantsPackage.from_IDs(['methane', 'ethane', 'propane', 'pentane',\
        'heptane', 'decane'])
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        #import pdb; pdb.set_trace()
        #constants.Vcs = Vc.tolist()
        constants.omegas = w.tolist()
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        kijs = Bin.tolist()

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)
        zs = [0.8097, 0.0566, 0.0306, 0.0457, 0.0330, 0.0244] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()
        #for i in range(100):
            #PT = flasher.flash(T=T[0], P=P[i], zs=zs)
        PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo: ' +str(dt)+' s')

        import pdb; pdb.set_trace()

    @unittest.skip("ok")
    def test_MY10(self):
        R = 8.3144598
        # C1, C2, C3, n-C4, n-C5, n-C6, n-C7, n-C8, n-C10, n-C14
        z = np.array([0.35, 0.03, 0.04, 0.06, 0.04, 0.03, 0.05, 0.05, 0.3, 0.05])[:,np.newaxis]*np.ones([10,100])
        Tc = np.array([190.6, 305.4, 369.8, 425.2, 469.6, 507.5, 540.3, 568.8, 617.9, 691.9]) # Kelvin
        Pc = np.array([45.4, 48.2, 41.9, 37.5, 33.3, 30.1, 27.4, 24.9, 21.0, 15.2])*101325 # pascal
        P = np.array([32.70])*101325*np.ones([100,])
        Pv = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        T = np.array([563.50])

        Mw = np.array([16.04, 30.07, 44.1, 58.12, 72.15, 86.178, 100.205, 114.232, 142.29, 198.39])*1e-3
        # Mw = np.array([16.04, 30.07, 44.1, 58.12, 72.15, 84, 107, 107, 147, 190])

        w = np.array([0.008, 0.098, 0.152, 0.193, 0.251, 0.305, 0.305, 0.396, 0.484, 0.747])

        #Bin = np.array([[0., 0.03], [0.03, 0.]])
        Bin = np.zeros([len(z), len(z)])

        Bin[0][3] = 0.02
        Bin[0][4] = 0.02
        Bin[0][5] = 0.025
        Bin[0][6] = 0.025
        Bin[0][7] = 0.035
        Bin[0][8] = 0.045
        Bin[0][9] = 0.045
        Bin[3][0] = 0.02
        Bin[4][0] = 0.02
        Bin[5][0] = 0.025
        Bin[6][0] = 0.025
        Bin[7][0] = 0.035
        Bin[8][0] = 0.045
        Bin[9][0] = 0.045

        CP1 = np.zeros_like(Tc)
        CP2 = np.zeros_like(Tc)
        CP3 = np.zeros_like(Tc)
        CP4 = np.zeros_like(Tc)

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)

        obj.run(z,Mw)

        ''' Using Thermo '''
        from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVL
        from thermo.interaction_parameters import IPDB
        import time
        # C1, C2, C3, n-C4, n-C5, n-C6, n-C7, n-C8, n-C10, n-C14
        constants, properties = ChemicalConstantsPackage.from_IDs(['methane', 'ethane', 'propane', 'butane',\
        'pentane', 'hexane', 'heptane', 'octane', 'decane', 'tetradecane'])
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        #import pdb; pdb.set_trace()
        #constants.Vcs = Vc.tolist()
        constants.omegas = w.tolist()
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        kijs = Bin.tolist()

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)
        zs = [0.35, 0.03, 0.04, 0.06, 0.04, 0.03, 0.05, 0.05, 0.3, 0.05] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()
        for i in range(100):
            PT = flasher.flash(T=T[0], P=P[i], zs=zs)
        #PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo: ' +str(dt)+' s')

        import pdb; pdb.set_trace()

    @unittest.skip("ok")
    def test_Dandekar(self):

        R = 10.73159
        z = np.array([0.8232,0.0871,0.0505,0.0198,0.0194])[:,np.newaxis]
        # C1, C3, n-C5, n-C10, n-C16
        Tc = np.array([343,666,845,1112,1291])#Rankine
        Pc = np.array([667.2,615.8,489.4,305.7,205.7]) #psia
        Bin = np.array([[0,0.009,0.021,0.052,0.080],\
                     [0.009,0,0.003,0.019,0.039],\
                     [0.021,0.003,0,0.008,0.022],\
                     [0.052,0.019,0.008,0,0.004],\
                     [0.08,0.039,0.022,0.004,0]])
        w = np.array([0.008,0.152,0.251,0.490,0.742])
        Mw = np.array([16.043,44.097,72.15,142.29,226.41])
        T = np.array([100+459.67]) #Fahrenheit
        P = np.array([1500]) #psia

        Pb_guess = 8e6

        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P,Pb_guess)
        obj.run(z,Mw)
        import pdb; pdb.set_trace()

    @unittest.skip("ok - poucos dados fornecidos, resultado aproximado")
    def test_Varavei_3_1_biphase(self):
        R = 8.3144598
        z = np.array([0.4, 0.1, 0.2, 0.3])[:,np.newaxis]*np.ones([4,100])
                     # H2O, C6, C10, C15
        T = np.array([333.15]) # Kelvin
        P = np.array([500])*6895*np.ones([100,]) # Pascal

        Tc = np.array([647.3, 507.4, 594.906, 676.266]) # Kelvin
        Pc = np.array([22089.00, 2969.00, 2439.00, 1824.00])*1000 # Pascal

        Bin = np.zeros([len(z), len(z)])

        w = np.array([0.344, 0.296, 0.5764, 0.7678])
        Mw = np.array([18, 86.2, 142.3, 206])*1e-3

        CP1 = np.zeros_like(Tc)
        CP2 = np.zeros_like(Tc)
        CP3 = np.zeros_like(Tc)
        CP4 = np.zeros_like(Tc)

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)
        obj.run(z,Mw)

        ''' Using Thermo '''
        from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVLN
        from thermo.interaction_parameters import IPDB
        import time
        # H2O, C6, C10, C15
        constants, properties = ChemicalConstantsPackage.from_IDs(['water', 'hexane', 'decane', 'pentadecane'])
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        #import pdb; pdb.set_trace()
        #constants.Vcs = Vc.tolist()
        constants.omegas = w.tolist()
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        kijs = Bin.tolist()

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVLN(constants, properties, liquids=[liquid, liquid], gas=gas)
        zs = [0.4, 0.1, 0.2, 0.3] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()
        for i in range(100):
            PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        #PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo: ' +str(dt)+' s')

        import pdb; pdb.set_trace()


    @unittest.skip("not ok - erro na simulação - Z negativos (dados suspeitos)")
    def test_Connolly_3_3_1(self):
        R = 8.3144598
                      # H20, C2-C11,   C12-C16,   C17-C21,   C22-C27,    C28-C35,   C36-C49,    C50+
        z = np.array([0.983415, 0.000325, 0.003012, 0.002570, 0.002252, 0.001912, 0.001676, 0.004838])[:,np.newaxis]
        Tc = np.array([647.37, 635.64, 701.24, 772.05, 826.30, 879.55, 936.97, 1260.0]) # Kelvin
        Pc = np.array([221.2, 24.11, 19.25, 15.10, 12.29, 9.94, 7.79, 6.00])*101325 # pascal
        P = np.array([85])*101325
        Pv = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        T = np.array([571])
        Mw = np.array([18.015, 143.66, 193.82, 263.40, 336.29, 430.48, 573.05, 1033.96])*1e-3
        w = np.array([0.344, 0.4645, 0.6087, 0.788, 0.9467, 1.1042, 1.273, 1.65])

        Bin = np.array([[0, 0.0952, -0.48, -0.48, 0.45, 0.53, 0.5, 0.5], \
                        [0.0952, 0, 0, 0, 0.13, 0.135, 0.1277, 0.1], \
                        [-0.48, 0, 0, 0, 0.05, 0.08, 0.1002, 0.1], \
                        [-0.48, 0, 0, 0, 0, 0, 0.09281, 0.130663], \
                        [0.45, 0.13, 0.05, 0, 0, 0, 0, 0.006], \
                        [0.53, 0.135, 0.08, 0, 0, 0, 0, 0.006], \
                        [0.5, 0.1277, 0.1002, 0.09281, 0, 0, 0, 0], \
                        [0.5, 0.1, 0.1, 0.130663, 0.006, 0.006, 0, 0]])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess)
        obj.run(z,Mw)

        import pdb; pdb.set_trace()


    @unittest.skip("ok")
    def test_Connolly332(self):
        R = 8.3144598
                    # H20,   C8,     C13,    C24,    C61+
        #z = np.array([0.5, 0.2227, 0.1402, 0.1016, 0.0355])[:,np.newaxis]
        z = np.array([0.5, 0.2227, 0.1402, 0.1016, 0.0355])[:,np.newaxis]*np.ones([5,1])
        Tc = np.array([647.37, 575.78, 698, 821.3, 1010.056]) # Kelvin
        Pc = np.array([221.2, 34.82, 23.37, 12.07, 7.79])*101325 # pascal

        #P = np.array([4.5])*101325
        P = np.array([20*101325])*np.ones([1,]) # Pascal

        Pv = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
        T = np.array([550])
        Mw = np.array([18.015, 116, 183, 337, 858])*1e-3
        w = np.array([0.344, 0.4, 0.84, 1.07, 1.33])

        Bin = np.array([[0, 0.5, 0.5, 0.5, 0.5], \
                        [0.5, 0, 0, 0, 0], [0.5, 0, 0, 0, 0], \
                        [0.5, 0, 0, 0, 0], [0.5, 0, 0, 0, 0]])

        CP1 = np.array([32.20, -1.23e+01, -5.08, -5.69, 0.1230])
        CP2 = np.array([0.001924, 6.65e-01, 9.97e-01, 1.840, 4.750])
        CP3 = np.array([1.055E-05, -2.52e-04, -4.14e-04, -7.64e-04, -1.95e-03])
        CP4 = np.array([-3.596e-09, 0.0, 0.0, 0.0, 0.0])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)
        obj.run(z,Mw)

        ''' Using Thermo '''
        from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVLN
        from thermo.interaction_parameters import IPDB
        import time
        # H20,   C8,     C13,    C24,    C61+
        constants, properties = ChemicalConstantsPackage.from_IDs(['water', 'octane', 'C13', 'C24', \
                                                                   'C60'])
        #import pdb; pdb.set_trace()
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        #import pdb; pdb.set_trace()
        #constants.Vcs = Vc.tolist()
        constants.omegas = w.tolist()
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        kijs = Bin.tolist()

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVLN(constants, properties, liquids=[liquid, liquid], gas=gas)
        zs = [0.5, 0.2227, 0.1402, 0.1016, 0.0355] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()
        for i in range(1):
            PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        #PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo: ' +str(dt)+' s')

        import pdb; pdb.set_trace()


    @unittest.skip("aparentemente ok - diagrama de fase")
    def test_Connolly333(self):
        R = 8.3144598
                    # H20, C6, C10, C15
        z = np.array([0.95, 0.005, 0.015, 0.03])[:,np.newaxis]
        Tc = np.array([647.3, 507.5, 622.1, 718.6]) # Kelvin
        Pc = np.array([220.4732, 32.88847, 25.34013, 18.49090])*101325 # pascal
        P = np.array([40])*101325
        Pv = np.array([0.0, 0.0, 0.0, 0.0])
        T = np.array([350])
        Mw = np.array([18.015, 86.178, 134.00, 206.00])*1e-3
        w = np.array([0.344, 0.275040, 0.443774, 0.651235])

        Bin = np.array([[0, 0.48, 0.48, 0.48], \
                        [0.48, 0, 0.002866, 0.010970], [0.48, 0.002866, 0, 0.002657], \
                        [0.48, 0.010970, 0.002657, 0]])

        CP1 = np.array([32.243, -4.413, -7.913, -11.916])
        CP2 = np.array([1.924e-3, 0.5820, 0.9609, 1.433])
        CP3 = np.array([1.055e-5, -3.119e-4, -5.288e-4, -7.972e-4])
        CP4 = np.array([-3.596e-9, 6.494e-8, 1.131e-7, 1.720e-7])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)
        obj.run(z,Mw)

        ''' Using Thermo '''
        from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVLN
        from thermo.interaction_parameters import IPDB
        import time
        # H20, C6, C10, C15
        constants, properties = ChemicalConstantsPackage.from_IDs(['water', 'C6', 'C10', 'C15'])
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        #import pdb; pdb.set_trace()
        #constants.Vcs = Vc.tolist()
        constants.omegas = w.tolist()
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        kijs = Bin.tolist()

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVLN(constants, properties, liquids=[liquid, liquid], gas=gas)
        zs = [0.95, 0.005, 0.015, 0.03] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()
        for i in range(1):
            PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        #PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo: ' +str(dt)+' s')

        import pdb; pdb.set_trace()


    @unittest.skip("ok - Diagrama de fase")
    def test_Connolly431(self):
        R = 8.3144598
        # H2O, PC1, PC2, PC3, PC4
        z = np.array([0.5, 0.15, 0.1, 0.1, 0.15])[:,np.newaxis]*np.ones([5,1])
        Tc = np.array([647.3, 305.556, 638.889, 788.889, 838.889]) # Kelvin
        Pc = np.array([220.8900, 48.82, 19.65, 10.20, 7.72])*101325 # pascal
        P = np.array([100*101325])*np.ones([1,])
        Pv = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        T = np.array([530])

        Mw = np.array([18.015, 30.00, 156.00, 310.00, 400.00])*1e-3
        w = np.array([0.344, 0.098, 0.535, 0.891, 1.085])

        Bin = np.array([[0, 0.71918, 0.45996, 0.26773, 0.24166], \
                        [0.71918, 0, 0, 0, 0], [0.45996, 0, 0, 0, 0], \
                        [0.26773, 0, 0, 0, 0], [0.24166, 0, 0, 0, 0]])

        CP1 = np.zeros_like(Tc)
        CP2 = np.zeros_like(Tc)
        CP3 = np.zeros_like(Tc)
        CP4 = np.zeros_like(Tc)

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)
        obj.run(z,Mw)

        import pdb; pdb.set_trace()

    @unittest.skip("aparentemente ok p/ o chines tb - não tem como certificar")
    def test_li2019_water_C3_C16(self):
        R = 8.3144598
                      # H20, C3,   C16
        z = np.array([0.75, 0.15, 0.1])[:,np.newaxis]
        Tc = np.array([647.3, 369.8, 717]) # Kelvin
        Pc = np.array([220.89, 42.46, 14.19])*101325 # pascal
        P = np.array([40])*101325
        Pv = np.array([0.0, 0.0, 0.0])
        T = np.array([460])
        Mw = np.array([18.02, 36.0321, 226.4412])*1e-3
        w = np.array([0.344, 0.152, 0.742])

        Bin = np.array([[0, 0.6841, 0.3583], [0.6841, 0, 0], [0.3583, 0, 0]])

        CP1 = np.array([32.243, -4.224, ])
        CP2 = np.array([1.924e-3, 0.3063, ])
        CP3 = np.array([1.055e-5, -1.586e-4, ])
        CP4 = np.array([-3.596e-9, 3.215e-8, ])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess)
        obj.run(z,Mw)

        ''' Using Thermo '''
        from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVLN
        from thermo.interaction_parameters import IPDB
        import time
        # H20, C3,   C16
        constants, properties = ChemicalConstantsPackage.from_IDs(['water', 'C3', 'C16'])
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        #import pdb; pdb.set_trace()
        #constants.Vcs = Vc.tolist()
        constants.omegas = w.tolist()
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        kijs = Bin.tolist()

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVLN(constants, properties, liquids=[liquid, liquid], gas=gas)
        zs = [0.75, 0.15, 0.1] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()
        for i in range(1):
            PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        #PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo: ' +str(dt)+' s')

        import pdb; pdb.set_trace()

    @unittest.skip("aparentemente ok p/ o chines tb - não tem como certificar")
    def test_li2019_water_C4_C20(self):
        R = 8.3144598
                      # H20, C4,  C20
        z = np.array([0.8, 0.16, 0.04])[:,np.newaxis]
        Tc = np.array([647.0, 425.2, 782.0]) # Kelvin
        Pc = np.array([220.5, 38.0, 14.6])*101325 # pascal
        P = np.array([60])*101325
        Pv = np.array([0.0, 0.0, 0.0])
        T = np.array([400])
        Mw = np.array([18.02, 58.12, 240.2140])*1e-3
        w = np.array([0.344, 0.1928, 0.8160])

        Bin = np.array([[0, 0.5, 0.5], [0.5, 0, 0], [0.5, 0, 0]])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess)
        obj.run(z,Mw)

        import pdb; pdb.set_trace()


    @unittest.skip("ok - Água/Vapor. ok para o chines")
    def test_Sofyan_case7(self):
        R = 8.3144598
                      # H2O / C1 / CO2 / H2S
        #z = np.array([0.1488, 0.2991, 0.0494, 0.5027])[:,np.newaxis]
        #z = np.array([0.1496, 0.3009, 0.0498, 0.4997])[:,np.newaxis]
        #z = np.array([0.0496, 0.0494, 0.4, 0.5])[:,np.newaxis]
        z = np.array([0.5008, 0.0504, 0.0503, 0.3986])[:,np.newaxis]
        Tc = np.array([647.3, 190.6, 304.2, 373.2]) # Kelvin
        Pc = np.array([220.5, 46, 73.8, 89.4])*101325 # pascal
        P = np.array([62.6])*101325
        Pv = np.array([0.0, 0.0, 0.0, 0.0])
        T = np.array([310.95])
        Mw = np.array([18.01528, 16.04, 44.01, 34.1])*1e-3
        w = np.array([0.344, 0.008, 0.225, 0.1])

        Bin = np.array([[0, 0.4928,  0, 0.04], \
                        [0.4928, 0,  0, 0], \
                        [0,  0, 0, 0], \
                        [0.04, 0,  0, 0]])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess)
        obj.run(z,Mw)

        import pdb; pdb.set_trace()

    @unittest.skip("ok - Água/Vapor")
    def test_Sabet_case1(self):
        R = 8.3144598
                    # H2O/ C1 / nC5 / nC10 / CO2 / H2S
        z = np.array([0.3, 0.45, 0.1, 0.05, 0.05, 0.05])[:,np.newaxis]
        Tc = np.array([647.3, 190.6, 469.6, 617.9, 304.2, 373.2]) # Kelvin
        Pc = np.array([220.5, 46, 33.3, 21, 73.8, 89.4])*101325 # pascal
        P = np.array([18e6])
        Pv = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        T = np.array([473.15])
        Mw = np.array([18.01528, 16.04, 60.05, 120.107, 44.01, 34.1])*1e-3
        w = np.array([0.344, 0.008, 0.251, 0.484, 0.225, 0.1])

        Bin = np.array([[0, 0.4907, 0.5, 0.45, 0.2, 0.275], \
                        [0.4907, 0, 0.0206, 0.0522, 0.103, 0.031], \
                        [0.5, 0.0206, 0, 0.0078, 0.125, 0.095], \
                        [0.45, 0.0522, 0.0078, 0, 0.11, 0.1], \
                        [0.2, 0.103, 0.125, 0.11, 0, 0.096], \
                        [0.275, 0.031, 0.095, 0.1, 0.096, 0]])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess)
        obj.run(z,Mw)

        import pdb; pdb.set_trace()

    @unittest.skip("ok - Água/Óleo")
    def test_Sabet_case2(self):
        R = 8.3144598
                    # H2O/ CO2 / H2S
        z = np.array([0.5, 0.45, 0.05])[:,np.newaxis]
        Tc = np.array([647.3, 304.2, 373.2]) # Kelvin
        Pc = np.array([220.5, 73.8, 89.4])*101325 # pascal
        P = np.array([15e6])
        Pv = np.array([0.0, 0.0, 0.0])
        T = np.array([333.15])
        Mw = np.array([18.01528, 44.01, 34.1])*1e-3
        w = np.array([0.344, 0.225, 0.1])

        Bin = np.array([[0, 0.2, 0.275], [0.2, 0, 0.096], [0.275, 0.096, 0]])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess)
        obj.run(z,Mw)

        import pdb; pdb.set_trace()

    @unittest.skip("ok - diagrama com 3 fases")
    def test_Sabet_case3(self):
        R = 8.3144598
                    # H2O / nC5 / nC10 / CO2 / H2S
        z = np.array([0.31, 0.03, 0.09, 0.42, 0.15])[:,np.newaxis]
        Tc = np.array([647.3, 469.6, 617.9, 304.2, 373.2]) # Kelvin
        Pc = np.array([220.5, 33.3, 21, 73.8, 89.4])*101325 # pascal
        P = np.array([13.1e6])
        Pv = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
        T = np.array([373.15])
        Mw = np.array([18.01528, 60.05, 120.107, 44.01, 34.1])*1e-3
        w = np.array([0.344, 0.251, 0.484, 0.225, 0.1])

        Bin = np.array([[0, 0.5, 0.45, 0.2, 0.275], \
                        [0.5, 0, 0.0078, 0.125, 0.095], \
                        [0.45, 0.0078, 0, 0.11, 0.1], \
                        [0.2, 0.125, 0.11, 0, 0.096], \
                        [0.275, 0.095, 0.1, 0.096, 0]])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess)
        obj.run(z,Mw)

        import pdb; pdb.set_trace()

    # 3 phases

    @unittest.skip("ok para a instabilidade em y")
    def test_Varavei_3_1_triphase(self):
        R = 8.3144598
                    # H2O, C1, C6, C10, C15
        #z = np.array([0.2, 0.1, 0.1, 0.2, 0.4])[:,np.newaxis]
        z = np.array([0.2, 0.1, 0.1, 0.2, 0.4])[:,np.newaxis]*np.ones([5,1])

        T = np.array([366.483]) # Kelvin
        P = np.array([200*6895.0])*np.ones([1,]) # Pascal
        #P = np.array([200*6895.0])

        Tc = np.array([647.3, 190.6, 507.4, 594.906, 676.266]) # Kelvin
        Pc = np.array([22089.00, 4600, 2969.00, 2439.00, 1824.00])*1000 # Pascal

        Bin = np.array([[0, 0, 0, 0, 0], \
                        [0, 0, 0, 0, 0], \
                        [0, 0, 0, 0, 0], \
                        [0, 0, 0, 0, 0], \
                        [0, 0, 0, 0, 0]])

        w = np.array([0.344, 0.008, 0.296, 0.5764, 0.7678])
        Mw = np.array([18, 16.04, 86.2, 142.3, 206])*1e-3

        # Dados do Stars
        CP1 = np.array([32.243, 19.251, -4.413, -7.913, -11.916])
        CP2 = np.array([1.924e-3, 5.213e-2, 0.5820, 0.9609, 1.433])
        CP3 = np.array([1.055e-5, 1.197e-5, -3.119e-4, -5.288e-4, -7.972e-4])
        CP4 = np.array([-3.596e-9, -1.132e-8, 6.494e-8, 1.131e-7, 1.720e-7])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)
        obj.run(z,Mw)

        ''' Using Thermo '''
        from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVLN
        from thermo.interaction_parameters import IPDB
        import time
        # H2O, C1, C6, C10, C15
        constants, properties = ChemicalConstantsPackage.from_IDs(['water', 'methane', 'hexane', 'decane', \
                                                                   'pentadecane'])
        #import pdb; pdb.set_trace()
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        #import pdb; pdb.set_trace()
        #constants.Vcs = Vc.tolist()
        constants.omegas = w.tolist()
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        kijs = Bin.tolist()

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVLN(constants, properties, liquids=[liquid, liquid], gas=gas)
        zs = [0.2, 0.1, 0.1, 0.2, 0.4] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()
        for i in range(1):
            PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        #PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo: ' +str(dt)+' s')

        import pdb; pdb.set_trace()

    @unittest.skip("ok vetorizado")
    def test_Sabet_case_3phases(self):
        R = 8.3144598
                    # H2O/ C1 / nC5 / nC10 / CO2 / H2S
        z = np.array([0.1, 0.3, 0.15, 0.25, 0.1, 0.1])[:,np.newaxis]
        z = np.array([0.1, 0.3, 0.15, 0.25, 0.1, 0.1])[:,np.newaxis]*np.ones([6,1000])

        Tc = np.array([647.3, 190.6, 469.6, 617.9, 304.2, 373.2]) # Kelvin
        Pc = np.array([220.5, 46, 33.3, 21, 73.8, 89.4])*101325 # pascal
        P = np.array([10e6])*np.ones([1000,])
        Pv = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        T = np.array([100 + 273.15])
        Mw = np.array([18.01528, 16.04, 60.05, 120.107, 44.01, 34.1])*1e-3
        w = np.array([0.344, 0.008, 0.251, 0.484, 0.225, 0.1])

        Bin = np.array([[0, 0.4907, 0.5, 0.45, 0.2, 0.275], \
                        [0.4907, 0, 0.0206, 0.0522, 0.103, 0.031], \
                        [0.5, 0.0206, 0, 0.0078, 0.125, 0.095], \
                        [0.45, 0.0522, 0.0078, 0, 0.11, 0.1], \
                        [0.2, 0.103, 0.125, 0.11, 0, 0.096], \
                        [0.275, 0.031, 0.095, 0.1, 0.096, 0]])

        CP1 = np.zeros_like(Tc)
        CP2 = np.zeros_like(Tc)
        CP3 = np.zeros_like(Tc)
        CP4 = np.zeros_like(Tc)

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)
        obj.run(z,Mw)

        ''' Using Thermo '''
        from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVLN
        from thermo.interaction_parameters import IPDB
        import time
        # H2O/ C1 / nC5 / nC10 / CO2 / H2S
        constants, properties = ChemicalConstantsPackage.from_IDs(['water', 'methane', 'pentane', 'decane', \
                                                                   'CO2', 'butane'])
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        #import pdb; pdb.set_trace()
        #constants.Vcs = Vc.tolist()
        constants.omegas = w.tolist()
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        kijs = Bin.tolist()

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVLN(constants, properties, liquids=[liquid, liquid], gas=gas)
        zs = [0.1, 0.3, 0.15, 0.25, 0.1, 0.1] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()
        for i in range(1000):
            PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        #PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo: ' +str(dt)+' s')

        import pdb; pdb.set_trace()

    @unittest.skip("ok - diagrama de fase - BUGANDO NO VETORIZADO PARA P=2 em T=400")
    def test_Connolly432(self):
        R = 8.3144598
                    # H2O,  CO2,      N2,       C1,      C2,      C3,     C4-C6,    PC1,      PC2,      PC3
        #z = np.array([0.1, 0.01089, 0.01746, 0.59391, 0.07821, 0.05319, 0.08703, 0.042705, 0.013635, 0.00297])[:,np.newaxis]
        z = np.array([0.1, 0.01089, 0.01746, 0.59391, 0.07821, 0.05319, 0.08703, 0.042705, 0.013635, 0.00297])[:,np.newaxis]*np.ones([10,1])
        Tc = np.array([647.3, 304.7, 126.2, 190.6, 305.43, 369.8, 448.08, 465.62, 587.8, 717.72]) # Kelvin
        Pc = np.array([220.8900, 73.86796, 33.94563, 46.04085, 48.83673, 42.65743, 35.50565, 28.32348, 17.06905, 11.06196])*101325 # pascal
        #P = np.array([100, 220, 80, 40, 200])*101325
        P = np.array([200])*101325
        Pv = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        #T = np.array([400])
        T = np.array([290.5])

        Mw = np.array([18.015, 44.01, 28.013, 16.043, 30.07, 44.097, 66.86942, 107.77943, 198.56203, 335.1979])*1e-3
        w = np.array([0.344, 0.225, 0.04, 0.013, 0.0986, 0.1524, 0.21575, 0.3123, 0.5567, 0.91692])

        Bin = np.array([[0, 0.0952, -0.48, -0.48, 0.45, 0.53, 0.50, 0.50, 0.50, 0.50], \
                        [0.0952, 0, 0, 0, 0.13, 0.135, 0.1277, 0.1, 0.1, 0.1], \
                        [-0.48, 0, 0, 0, 0.05, 0.08, 0.1002, 0.1, 0.1, 0.1], \
                        [-0.48, 0, 0, 0, 0, 0, 0.09281, 0.130663, 0.130663, 0.130663], \
                        [0.45, 0.13, 0.05, 0, 0, 0, 0, 0.006, 0.006, 0.006], \
                        [0.53, 0.135, 0.08, 0, 0, 0, 0, 0.006, 0.006, 0.006], \
                        [0.50, 0.1277, 0.1002, 0.09281, 0, 0, 0, 0, 0, 0], \
                        [0.50, 0.1, 0.1, 0.130663, 0.006, 0.006, 0, 0, 0, 0], \
                        [0.50, 0.1, 0.1, 0.130663, 0.006, 0.006, 0, 0, 0, 0], \
                        [0.50, 0.1, 0.1, 0.130663, 0.006, 0.006, 0, 0, 0, 0]])

        CP1 = np.zeros_like(Tc)
        CP2 = np.zeros_like(Tc)
        CP3 = np.zeros_like(Tc)
        CP4 = np.zeros_like(Tc)

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)

        obj.run(z,Mw)

        ''' Using Thermo '''
        from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVLN
        from thermo.interaction_parameters import IPDB
        import time
        # H2O,  CO2,      N2,       C1,      C2,      C3,     C4-C6,    PC1,      PC2,      PC3
        constants, properties = ChemicalConstantsPackage.from_IDs(['water', 'CO2', 'CO2', \
                                                                   'methane', 'ethane', 'propane', 'butane', \
                                                                   'ethane', 'propane', 'butane'])
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        #import pdb; pdb.set_trace()
        #constants.Vcs = Vc.tolist()
        constants.omegas = w.tolist()
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        kijs = Bin.tolist()

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVLN(constants, properties, liquids=[liquid, liquid], gas=gas)
        zs = [0.1, 0.01089, 0.01746, 0.59391, 0.07821, 0.05319, 0.08703, 0.042705, 0.013635, 0.00297] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()
        for i in range(1):
            PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        #PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo: ' +str(dt)+' s')

        import pdb; pdb.set_trace()

    @unittest.skip("mais ou menos ok - Vetorizado pau no flash bifasico")
    def test_lapene_water_benzene_toluene(self):
        R = 8.3144598
                      # H20 / tolueno / benzeno
        z = np.array([0.29, 0.01, 0.7])[:,np.newaxis]*np.ones([3,2])
        Tc = np.array([647, 593, 562]) # Kelvin
        Pc = np.array([220.5, 41, 48.9])*101325 # pascal
        P = np.array([1*101325])*np.ones([2,])
        Pv = np.array([0.0, 0.0, 0.0])
        T = np.array([345])
        Mw = np.array([18.02, 92.13, 78.11])*1e-3
        w = np.array([0.344, 0.262, 0.212])

        Bin = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess)
        obj.run(z,Mw)

        import pdb; pdb.set_trace()

    @unittest.skip("ok - diagrama de fase - estabilidade em P=50 nao detecta a terceira fase - pau no flash bifasico para o vetorizado")
    def test_lapene_water_nitrogen_C10_C20(self):
        R = 8.3144598
                      # H20 / nitrogen / C10 / C20
        z = np.array([0.55, 0.1, 0.1, 0.25])[:,np.newaxis]*np.ones([4,2])
        Tc = np.array([647, 126.2, 622, 782]) # Kelvin
        Pc = np.array([220.5, 34, 25.3, 14.6])*101325 # pascal
        #P = np.array([190])*101325
        P = np.array([40, 80])*101325
        Pv = np.array([0.0, 0.0, 0.0, 0.0])
        T = np.array([450])
        Mw = np.array([18, 28, 134, 275])*1e-3
        w = np.array([0.344, 0.04, 0.443, 0.816])

        Bin = np.array([[0, 0, 0, 0], [0, 0, 0, 0], \
                        [0, 0, 0, 0], [0, 0, 0, 0]])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess)
        obj.run(z,Mw)

        import pdb; pdb.set_trace()

    @unittest.skip("ok vetorizado")
    def test_lapene_18_components(self):
        R = 8.3144598

        z = np.array([20, 0.26, 3.60, 74.12, 7.94, 3.29, 0.68, 1.24, 0.55, 0.61, \
                    0.87, 1.15, 1.07, 0.95, 0.67, 1.65, 1.13, 0.22])[:,np.newaxis]
        z = z/np.sum(z)
        #z = z*np.ones([18,8])

        Tc = np.array([647.37, 126.20, 304.21, 190.60, 305.40, 369.80, 408.10, 425.20, \
                        464.74, 469.60, 515.28, 553.84, 581.28, 609.35, 626.97, 658.15, 778.15, 998.15]) # Kelvin
        Pc = np.array([221.20, 33.94, 73.77, 46.00, 48.84, 42.46, 36.48, 38.00, 34.77, \
                        33.74, 32.57, 31.00, 28.50, 26.50, 24.60, 21.20, 15.70, 13.50])*101325 # pascal
        #P = np.array([100, 100, 200, 200, 300, 300, 500, 500])*101325
        P = np.array([200])*101325

        Pv = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        T = np.array([450])
        Mw = np.array([18.02, 28.01, 44.01, 16.04, 30.07, 44.10, 58.12, 58.12, 71.94, \
                        72.15, 84.99, 97.87, 111.54, 126.10, 140.14, 179.30, 290.60, 450.00])*1e-3
        w = np.array([0.344, 0.04, 0.2250, 0.0115, 0.0908, 0.1454, 0.1760, 0.1928, 0.2235, 0.2273, \
                        0.2637, 0.2897, 0.3245, 0.3791, 0.4363, 0.5200, 0.6500, 0.7200])

        Bin = np.array([[0, 0.4778, 0.1896, 0.4850, 0.4920, 0.5525, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], \
                        [0.4778, 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1], \
                        [0.1896, 0.1, 0, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12], \
                        [0.485, 0.1, 0.12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.02, 0.03, 0.05, 0.07, 0.085, 0.07], \
                        [0.492, 0.1, 0.12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.02, 0.04, 0.05, 0.04], \
                        [0.5525, 0.1, 0.12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0.1, 0.12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0.1, 0.12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0.1, 0.12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0.1, 0.12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0.1, 0.12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0.1, 0.12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0.1, 0.12, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0.1, 0.12, 0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0.1, 0.12, 0.05, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0.1, 0.12, 0.07, 0.04, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0.1, 0.12, 0.085, 0.05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0.1, 0.12, 0.07, 0.04, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

        CP1 = np.zeros_like(Tc)
        CP2 = np.zeros_like(Tc)
        CP3 = np.zeros_like(Tc)
        CP4 = np.zeros_like(Tc)

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)
        obj.run(z,Mw)

        ''' Using Thermo '''
        from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVLN
        from thermo.interaction_parameters import IPDB
        import time

        constants, properties = ChemicalConstantsPackage.from_IDs(['water', 'methane', 'hexane', 'decane', \
                                                                   'water', 'methane', 'hexane', 'decane', \
                                                                   'water', 'methane', 'hexane', 'decane', \
                                                                   'water', 'methane', 'hexane', 'decane', \
                                                                   'water', 'methane'])
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        #import pdb; pdb.set_trace()
        #constants.Vcs = Vc.tolist()
        constants.omegas = w.tolist()
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        kijs = Bin.tolist()

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVLN(constants, properties, liquids=[liquid, liquid], gas=gas)
        zs = [0.16666667, 0.00216667, 0.03, 0.61766667, 0.06616667, 0.02741667, 0.00566667, \
              0.01033333, 0.00458333, 0.00508333, 0.00725, 0.00958333, 0.00891667, 0.00791667, \
              0.00558333, 0.01375, 0.00941667, 0.00183333] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()
        for i in range(100):
            PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        #PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo: ' +str(dt)+' s')


        import pdb; pdb.set_trace()

    @unittest.skip("not ok - região trifásica muito estreita, pau no flash bifasico")
    def test_lapene_4_5(self):
        R = 8.3144598

        z = np.array([0.983415, 0.000325, 0.003012, 0.002570, 0.002252, 0.001912, 0.001676, 0.004838])[:,np.newaxis]
        Tc = np.array([647.37, 635.64, 701.24, 772.05, 826.30, 879.55, 936.97, 1260.00]) # Kelvin
        Pc = np.array([221.20, 24.11, 19.25, 15.10, 12.29, 9.94, 7.79, 6.00])*101325 # pascal
        P = np.array([7])*101325
        Pv = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        T = np.array([438.5])
        Mw = np.array([18.02, 143.66, 193.82, 263.40, 336.29, 430.48, 573.05, 1033.96])*1e-3
        w = np.array([0.344, 0.4645, 0.6087, 0.7880, 0.9467, 1.1042, 1.2730, 1.6500])

        Bin = np.array([[0, 0, 0, 0, 0, 0, 0, 0], \
                        [0, 0, 0, 0, 0, 0, 0, 0], \
                        [0, 0, 0, 0, 0, 0, 0, 0], \
                        [0, 0, 0, 0, 0, 0, 0, 0], \
                        [0, 0, 0, 0, 0, 0, 0, 0], \
                        [0, 0, 0, 0, 0, 0, 0, 0], \
                        [0, 0, 0, 0, 0, 0, 0, 0], \
                        [0, 0, 0, 0, 0, 0, 0, 0]])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess)
        obj.run(z,Mw)

        import pdb; pdb.set_trace()

    @unittest.skip("not ok")
    def test_heidari_water_C6_C10_C15(self):
        R = 8.3144598
                      # H20 / C6 / C10 / C15
        z = np.array([0.55, 0.1, 0.1, 0.25])[:,np.newaxis]
        Tc = np.array([647.3, 507.4, 594.906, 676.266]) # Kelvin
        Pc = np.array([22089.00, 2969.00, 2439.00, 1824.00])*1000 # pascal
        P = np.array([500])*1000
        Pv = np.array([0.0, 0.0, 0.0, 0.0])
        T = np.array([388])
        Mw = np.array([18, 86.178, 142.3, 206])*1e-3
        w = np.array([0.344, 0.296, 0.5764, 0.7678])

        Bin = np.array([[0, 0.48, 0.48, 0.48], [0.48, 0, 0.00280, 0.01097], \
                        [0.48, 0.00280, 0, 0.02657], [0.48, 0.01097, 0.02657, 0]])

        CP1 = np.zeros_like(Tc)
        CP2 = np.zeros_like(Tc)
        CP3 = np.zeros_like(Tc)
        CP4 = np.zeros_like(Tc)

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)

        obj.run(z,Mw)

        import pdb; pdb.set_trace()

    @unittest.skip("not ok")
    def test_heidari_case_4_7_1_1(self):
        R = 8.3144598
                      # H20 / N2 / C10 / C20
        z = np.array([0.55, 0.1, 0.1, 0.25])[:,np.newaxis]
        Tc = np.array([647.3, 126.2, 622.0, 782.0]) # Kelvin
        Pc = np.array([22050.0, 3400.0, 2530.0, 1460.0])*1000 # pascal
        P = np.array([50])*101325
        Pv = np.array([0.0, 0.0, 0.0, 0.0])
        T = np.array([450])
        Mw = np.array([18, 28.0134, 142.3, 240.2140])*1e-3
        w = np.array([0.344, 0.04, 0.433, 0.816])

        Bin = np.array([[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], \
                        [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]])

        CP1 = np.zeros_like(Tc)
        CP2 = np.zeros_like(Tc)
        CP3 = np.zeros_like(Tc)
        CP4 = np.zeros_like(Tc)

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)

        obj.run(z,Mw)

        import pdb; pdb.set_trace()

    @unittest.skip("ok - diagrama de fase")
    def test_heidari_water_C10_C15_C19_C25_C26(self):
        R = 8.3144598
                      # H20 / C10 / C15 / C19 / C25 / C26
        z = np.array([0.0915, 0.0004, 0.0009, 0.0009, 0.0023, 0.9040])[:,np.newaxis]
        Tc = np.array([647.3, 594.906, 676.266, 726.653, 791.757, 802.131]) # Kelvin
        Pc = np.array([22089.00, 2439.00, 1824.00, 1581.00, 1435.00, 1420.00])*1000 # pascal
        P = np.array([1500])*1000
        Pv = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        T = np.array([500])
        Mw = np.array([18, 142.3, 206, 268.31, 300.2675, 412.7])*1e-3
        w = np.array([0.344, 0.5764, 0.7678, 0.9046, 1.0755, 1.1014])

        Bin = np.array([[0, 0.48, 0.48, 0.48, 0.48, 0.48], [0.48, 0, 0, 0, 0, 0], \
                        [0.48, 0, 0, 0, 0, 0], [0.48, 0, 0, 0, 0, 0], \
                        [0.48, 0, 0, 0, 0, 0], [0.48, 0, 0, 0, 0, 0]])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess)
        obj.run(z,Mw)

        import pdb; pdb.set_trace()


    @unittest.skip("ok - diagrama de fase")
    def test_heidari_4_8_1_3(self):
        R = 8.3144598
                      # H20 / CO2 / C2 / nC5 / C8 / C12
        z = np.array([0.134, 0.800, 0.033, 0.011, 0.011, 0.011])[:,np.newaxis]
        Tc = np.array([647.3, 304.2, 305.4, 469.6, 570.5, 663.9]) # Kelvin
        Pc = np.array([22048.3, 7376.1, 4883.6, 3374.0, 2950.4, 2191.6])*1000 # pascal
        P = np.array([7040])*1000
        Pv = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        T = np.array([460])
        Mw = np.array([18, 44.01, 30.070, 60.05, 96.08, 162.27])*1e-3
        w = np.array([0.344, 0.2250, 0.098, 0.251, 0.351, 0.522])

        Bin = np.array([[0, 0.2, 0.4911, 0.5, 0.48, 0.48], \
                        [0.2, 0, 0.15, 0.15, 0.15, 0.15], \
                        [0.4911, 0.15, 0, 0.008578, 0.01796, 0.033751], \
                        [0.5, 0.15, 0.008578, 0, 0.001765, 0.008637], \
                        [0.48, 0.15, 0.01796, 0.001765, 0, 0.002618], \
                        [0.48, 0.15, 0.033751, 0.008637, 0.002618, 0]])

        CP1 = np.array([33.75536, 29.26153, 33.30586, 33.77337, -20.53700, -26.45610])
        CP2 = np.array([-0.00594, -0.02236, -0.01113, 0.24845, 0.69221, 1.02055])
        CP3 = np.array([2.24e-05, 0.000265, 0.000357, 0.000253, -0.000280, -0.000410])
        CP4 = np.array([-9.96e-09, -4.15e-07, -3.76e-07, -3.84e-07, 0, 0])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)
        obj.run(z,Mw)

        ''' Using Thermo '''
        from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVLN
        from thermo.interaction_parameters import IPDB
        import time
        # H20 / CO2 / C2 / nC5 / C8 / C12
        constants, properties = ChemicalConstantsPackage.from_IDs(['water', 'CO2', 'ethane', 'pentane', 'octane', 'dodecane'])
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        #constants.Vcs = Vc.tolist()
        constants.omegas = w.tolist()
        #kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        kijs = Bin.tolist()
        #kijs = Bin

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVLN(constants, properties, liquids=[liquid, liquid], gas=gas)
        zs = [0.134, 0.800, 0.033, 0.011, 0.011, 0.011] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()
        for i in range(1):
            #print(i)
            PT = flasher.flash(T=T[i], P=P[i], zs=z.T[i])
        #PT = flasher.flash(T=T, P=P, zs=z)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo: ' +str(dt)+' s')

        import pdb; pdb.set_trace()


    @unittest.skip("ok - diagrama de fase")
    def test_heidari_4_8_1_5(self):
        R = 8.3144598
                      # H20 / C1 / C40
        z = np.array([0.3333, 0.3333, 0.3334])[:,np.newaxis]
        Tc = np.array([647.3, 190.6, 934.3]) # Kelvin
        Pc = np.array([22048.32, 4599.93, 800.43])*1000 # pascal
        P = np.array([1000])*1000
        Pv = np.array([0.0, 0.0, 0.0])
        T = np.array([350])
        Mw = np.array([18, 12.0107, 480.4280])*1e-3
        w = np.array([0.344, 0.008, 1.259])

        Bin = np.array([[0, 0.4907, 0.4800], \
                        [0.4907, 0, 0.125465], \
                        [0.48, 0.125465, 0]])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess)
        obj.run(z,Mw)

        import pdb; pdb.set_trace()

    @unittest.skip("not ok")
    def test_Connolly_5_4_2(self):
        R = 8.3144598
                      #
        z = np.array([0.800, 0.066, 0.022, 0.022, 0.022, 0.068])[:,np.newaxis]
        Tc = np.array([304.20, 190.60, 369.80, 469.60, 568.80, 658.30]) # Kelvin
        Pc = np.array([73.76, 46.00, 42.46, 33.74, 24.82, 18.24])*101325 # pascal
        P = np.array([77.5])*101325
        Pv = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        T = np.array([303.35])
        Mw = np.array([44.01, 16.043, 44.097, 72.151, 114.232, 170.34])*1e-3
        w = np.array([0.225, 0.008, 0.152, 0.251, 0.394, 0.562])

        Bin = np.array([[0, 0.12, 0.12, 0.12, 0.1, 0.1], \
                        [0.12, 0, 0, 0, 0.0496, 0], \
                        [0.12, 0, 0, 0, 0, 0], \
                        [0.12, 0, 0, 0, 0, 0], \
                        [0.1, 0.0496, 0, 0, 0, 0], \
                        [0.1, 0, 0, 0, 0, 0]])

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess)
        obj.run(z,Mw)

        import pdb; pdb.set_trace()

    @unittest.skip("ok")
    def test_Aug_water_NWE(self):
        R = 8.3144598
                     # H2O, CO2,       C1,      C2_3,  C4_6,   C7_14,    C15_24,    C25
        z = np.array([0.5, 0.251925, 0.050625, 0.0295, 0.0371, 0.071575, 0.03725, 0.022025])[:,np.newaxis] # Mol %
        #z = z**np.ones([8,1]) # Número de componentes x número de bloco do pressão

        Tc = np.array([647.35, 304.20, 190.60, 343.64, 466.41, 603.07, 733.79, 923.20]) # Kelvin
        Pc = np.array([221.00, 73.76, 46.00, 45.05, 33.50, 24.24, 18.03, 17.26])*101325 # pascal
        P = np.array([400])*101325

        Pv = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        T = np.array([570])

        Mw = np.array([18.00, 44.01, 16.04, 37.2, 69.5, 140.96, 280.99, 519.62])*1e-3 # Paper Full EOS
        w = np.array([0.3434, 0.225, 0.008, 0.130, 0.244, 0.600, 0.903, 1.229])

        Bin = np.array([[0, 0.1896, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], \
                        [0.1896, 0, 0.12, 0.12, 0.12, 0.09, 0.09, 0.09], \
                        [0.4850, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0, 0, 0, 0, 0, 0, 0]])

        CP1 = np.zeros_like(Tc)
        CP2 = np.zeros_like(Tc)
        CP3 = np.zeros_like(Tc)
        CP4 = np.zeros_like(Tc)

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)

        obj.run(z,Mw)

        ''' Using Thermo '''
        from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVLN
        from thermo.interaction_parameters import IPDB
        import time
        # H2O, CO2,       C1,      C2_3,  C4_6,   C7_14,    C15_24,    C25
        constants, properties = ChemicalConstantsPackage.from_IDs(['water', 'CO2', 'methane', 'ethane', 'propane', \
                                                                   'hexadecane', 'hexadecane', 'hexadecane'])
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        #import pdb; pdb.set_trace()
        #constants.Vcs = Vc.tolist()
        constants.omegas = w.tolist()
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        kijs = Bin.tolist()

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVLN(constants, properties, liquids=[liquid, liquid], gas=gas)
        zs = [0.5, 0.251925, 0.050625, 0.0295, 0.0371, 0.071575, 0.03725, 0.022025] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()
        for i in range(1):
            PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        #PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo: ' +str(dt)+' s')

        import pdb; pdb.set_trace()

    @unittest.skip("not ok - dados de entrada confusos, teste de william")
    def test_Aug_water_BSB(self):
        R = 8.3144598

        z = np.array([0.75, 0.105055, 0.012915, 0.022545, 0.025065, 0.04956, 0.024165, 0.010695])[:,np.newaxis] # Mol %
        #z = z/np.sum(z)
        #z = z**np.ones([8,1]) # Número de componentes x número de bloco do pressão

        Tc = np.array([647.35, 304.22, 160.00, 344.22, 463.22, 605.78, 751.00, 942.50]) # Kelvin
        Pc = np.array([221.00, 73.77, 46.00, 45.00, 34.00, 21.75, 16.54, 16.42])*101325 # pascal
        P = np.array([200])*101325

        Pv = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        T = np.array([500])

        Mw = np.array([18.00, 44.01, 16.04, 37.2, 69.5, 140.96, 280.99, 519.62])*1e-3 # Paper Full EOS
        w = np.array([0.3434, 0.225, 0.008, 0.131, 0.240, 0.618, 0.957, 1.268])

        Bin = np.array([[0, 0.1896, 0.4850, 0.5, 0.5, 0.5, 0.5, 0.5], \
                        [0.1896, 0, 0.055, 0.055, 0.055, 0.105, 0.105, 0.105], \
                        [0.4850, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0, 0, 0, 0, 0, 0, 0], \
                        [0.5, 0, 0, 0, 0, 0, 0, 0]])

        CP1 = np.zeros_like(Tc)
        CP2 = np.zeros_like(Tc)
        CP3 = np.zeros_like(Tc)
        CP4 = np.zeros_like(Tc)

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)

        obj.run(z,Mw)

        import pdb; pdb.set_trace()

    @unittest.skip("not ok - testar")
    def test_Huang_Yang_case2(self):
        R = 8.3144598
                     # H2O, CO2, C16H34
        z = np.array([0.75, 0.15, 0.10])[:,np.newaxis] # Mol %
        #z = z**np.ones([8,1]) # Número de componentes x número de bloco do pressão

        Tc = np.array([647.3, 304.2, 717.0]) # Kelvin
        Pc = np.array([220.89, 73.76, 14.19])*101325 # pascal
        P = np.array([80])*101325

        Pv = np.array([0.0, 0.0, 0.0])
        T = np.array([500])

        Mw = np.array([18.0, 44.01, 226.41])*1e-3
        w = np.array([0.344, 0.225, 0.742])

        Bin = np.array([[0.0, 0.1896, 0.3583], \
                        [0.1896, 0.0, 0.125], \
                        [0.3583, 0.125, 0.0]])

        CP1 = np.array([32.20, 19.795, -13.00]) # J/(mol * Kelvin)
        CP2 = np.array([0.001907, 0.07343, 1.529]) # J/(mol * Kelvin**2)
        CP3 = np.array([1.06e-5, -5.602e-5, -8.5e-4]) # J/(mol * Kelvin**3)
        CP4 = np.array([-3.596e-9, 1.715e-8, 1.850e-7]) # J/(mol * Kelvin**4)

        Pb_guess = 8e6
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)

        obj.run(z,Mw)

        import pdb; pdb.set_trace()

    @unittest.skip("ok")
    def test_Aug_water_411(self):
        R = 8.3144598
                    # H2O, CO2, C1, C16
        z = np.array([0.3, 0.3, 0.1, 0.3])[:,np.newaxis]*np.ones([4,10000])
        #z = z**np.ones([8,1]) # Número de componentes x número de bloco do pressão

        Tc = np.array([647.3, 304.20, 190.60, 717.00]) # Kelvin
        Pc = np.array([220.48, 73.76, 45.96, 14.19])*101325 # pascal
        P = np.array([100])*101325*np.ones([10000,])

        Pv = np.array([0.0, 0.0, 0.0, 0.0])
        T = np.array([450])*np.ones([10000,])

        Mw = np.array([18.00, 44.01, 16.043, 226.4412])*1e-3
        w = np.array([0.344, 0.225, 0.008, 0.742])

        Bin = np.array([[0.0, 0.1896, 0.4850, 0.5000], \
                        [0.1896, 0.0, 0.1, 0.1250], \
                        [0.4850, 0.1, 0.0, 0.0780], \
                        [0.5, 0.1250, 0.0780, 0.0]])

        CP1 = np.zeros_like(Tc)
        CP2 = np.zeros_like(Tc)
        CP3 = np.zeros_like(Tc)
        CP4 = np.zeros_like(Tc)

        Pb_guess = 8e6
        #obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, Pb_guess, CP1, CP2, CP3, CP4)

        #obj.run(z,Mw)

        ''' Using Thermo '''
        from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVLN
        from thermo.interaction_parameters import IPDB
        import time
        # H2O, CO2, C1, C16
        constants, properties = ChemicalConstantsPackage.from_IDs(['water', 'CO2', 'methane', 'hexadecane'])
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        #constants.Vcs = Vc.tolist()
        constants.omegas = w.tolist()
        #kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        #kijs = Bin.tolist()
        kijs = Bin

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVLN(constants, properties, liquids=[liquid, liquid], gas=gas)
        zs = [0.3, 0.3, 0.1, 0.3] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()
        for i in range(10000):
            #print(i)
            PT = flasher.flash(T=T[i], P=P[i], zs=z.T[i])
        #PT = flasher.flash(T=T, P=P, zs=z)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo: ' +str(dt)+' s')

        import pdb; pdb.set_trace()

    

import numpy as np
from stability_check import StabilityCheck
import unittest
import matplotlib.pyplot as plt
from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVL
from thermo.interaction_parameters import IPDB
import time


class Li_et_al_table4:
    def __init__(self):
        self.R = 10.73159
        self.Tc = np.array([-116.59, 305.69, 453.65, 652.01])
        self.Tc = (self.Tc + 459.67)
        self.Pc = np.array([667.2,551.1,430.59,305.68])
        self.w = np.array([0.008,0.193,0.296,0.49])
        self.Mw = np.array([16.043,58.124,86.16,142.29])
        self.P = np.array([600]) * np.ones(5)
        self.T = np.array([160 + 459.67])
        self.Bin = np.zeros([len(self.w),len(self.w)])

class testes_Li_et_al_table4(unittest.TestCase):

    #Units: English system
    # methane n-butane n-hexane n-decane
    @unittest.skip("ok")
    def test_all(self):
        prop = Li_et_al_table4()
        z = np.array([[-0.58,0.38,0.6,0.6], [0,0,0,1], [0,0,1,0], [0,1,0,0], [1,0,0,0]])
        obj = StabilityCheck(prop.w,prop.Bin,prop.R,prop.Tc,prop.Pc,prop.T,prop.P)

        obj.run(z.T,prop.Mw)
        x = np.array([[0.1676,0.2239,0.3088,0.2997],[0.1724,0.0000,0.0000,0.8276], \
                    [0.1789,0.0000,0.8211,0.0000], [0.1639,0.8361,0.0000,0.0000],
                    [0.1724,0.0000,0.0000,0.8276]])
        y = np.array([[0.9115,0.0686,0.0190,0.0009], [0.9980,0.0000,0.0000,0.0020],
                    [0.9507,0.0000,0.0493,0.0000], [0.6960,0.3040,0.0000,0.0000],
                    [0.9980,0.0000,0.0000,0.0020]])

        print('x: ', 'y: ', obj.x, obj.y)
        print('L: ','V: ', obj.L, obj.V)

        for j in range(5):
            for i in range(obj.Nc):
                self.assertAlmostEqual(obj.x[i,j],x.T[i,j],4,'ValueError:Failed')
                self.assertAlmostEqual(obj.y[i,j],y.T[i,j],4,'ValueError:Failed')

class testes_casos(unittest.TestCase):
    @unittest.skip("ok")
    def teste_caso_Igor(self):
        R = 8.3144598
        z = np.array([0.35, 0.03, 0.04, 0.06, 0.04, 0.03, 0.05, 0.05, 0.3, 0.05])[:,np.newaxis]#*np.ones([10]) #exemplo aleatório
        P = np.array([104.9])*101325
        #P = np.array([14619094.696946101])*np.ones(10)
        #T = np.array([15 + 273.15])
        T = np.array([509.1])
        Tc = np.array([190.6, 305.4, 369.8, 425.2, 469.6 ])
        Pc =np.array([2109795.64])
        w = np.array([0.489])
        Bin = np.array([[0]])
        Mw = np.array([142.28e-3])
        C7 = np.array([0])#*np.ones([1,10])
        print('\ncaso5:')
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, 9e6)
        obj.run(z,Mw)

        print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L,'V: ',obj.V)
        #print('fl: ',obj.fl,'fv: ',obj.fv)

    @unittest.skip("ok")
    def teste_caso5(self):
        R = 8.3144598
        z = np.array([1.])[:,np.newaxis]#*np.ones([10]) #exemplo aleatório
        P = np.array([13.78951458E6])
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
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, 9e6)
        obj.run(z,Mw)

        print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L,'V: ',obj.V)
        #print('fl: ',obj.fl,'fv: ',obj.fv)

    @unittest.skip("?")
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

        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P, 9e6)
        obj.run(z,Mw)
        plt.plot(P/6894.757, -obj.dVtdP * 6894.757)
        plt.savefig('teste.png')
        import pdb; pdb.set_trace()

    @unittest.skip("?")
    def test_Firoozabadi(self):
        #methane, n-butane and n-decane
        R = 8.3144598
        z = np.array([0.9, 0.1, 0.065])[:,np.newaxis] * np.ones((3,100))
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
        z = np.array([0.95, 0.05, 0.0])[:,np.newaxis]*np.ones([3,100])
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
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P,Pv)
        obj.run(z,Mw)

        ''' Using Thermo '''
        constants, properties = ChemicalConstantsPackage.from_IDs(['methane', 'propane', 'hexane',\
        'decane', 'pentadecane', 'icosane'])
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        constants.Vcs = vc.tolist()
        constants.omegas = w.tolist()
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        kijs = Bin.tolist()

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVLN(constants, properties, liquid=liquid, gas=gas)
        zs = [0.95, 0.05, 0.0] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()
        PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo: ' +str(dt)+' s')

        t0_thermo = time.time()
        PT_vec = flasher.grid_flash(Ts=T, Ps=P, zs=z.T)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo vectorized: ' +str(dt)+' s')
        import pdb; pdb.set_trace()

    def teste_inj_fluid_ex1_firoo(self):
        R = 8.3144598
        z = np.array([1.0,0])[:,np.newaxis]
        Tc = np.array([190.56, 369.83])
        Pc = np.array([4599000, 4248000])
        vc = np.array([9.85988846e-05, 2.00001177e-04])
        #T = np.array([273.15])
        P = np.array([1e5])
        T = np.array([293])
        Mw = np.array([16.043e-3, 44.096e-3])
        w = np.array([0.0115, 0.1523])
        Bin = np.array([[0.0,0.],[0.,0.0]])
        print('\nExemplo firooz 2k:')

        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P)
        obj.run(z,Mw)
        import pdb; pdb.set_trace()

        constants, properties = ChemicalConstantsPackage.from_IDs(['methane', 'propane'])
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        constants.Vcs = vc.tolist()
        constants.omegas = w.tolist()
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        kijs = Bin.tolist()

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)
        zs = [1,0] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()
        PT = flasher.flash(T=T[0], P=P[0], zs=zs)
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

    @unittest.skip("ok")
    def test_water_inj_schmall3d(self):
        #  C1 C3 C6 C10 C15 C20
        R = 8.3144598
        z = np.array([0.5, 0.03, 0.07, 0.2, 0.15, 0.05])[:,np.newaxis]#* np.ones((6,100))
        Tc = np.array([190.6, 369.8, 507.4, 617.6, 708, 768])
        Pc = np.array([4600155, 4245517.5, 2968822.5, 2107560.0, 1.47e6, 1.17e6])
        vc = np.array([0.000099, 0.000203, 0.00037, 0.000603, 0.000895, 0.00169])
        Mw = np.array([16.042e-3, 44.1e-3, 86.178e-3, 142.276e-3, 212.41e-3, 282.5e-3])
        w = np.array([0.008, 0.152, 0.299, 0.489, 0.685, 0.912])
        Bin = np.array([[0.,0.,.0,0.,0.,0.], [0.,0.,.0,0.,0.,0.], [0.,0.,.0,0.,0.,0.], [0.,0.,.0,0.,0.,0.],
                        [0.,0.,.0,0.,0.,0.], [0.,0.,.0,0.,0.,0.]])
        Pv = np.array([8e6, 0., 0., 0., 0., 0.])
        #P = np.linspace(8.46e6, 62e6, 100)
        P = np.array([10.34e6])#*np.ones(100)
        T = np.array([344.25])#*np.ones(100)

        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P)
        obj.run(z,Mw)

        ''' Using Thermo '''
        constants, properties = ChemicalConstantsPackage.from_IDs(['methane', 'propane', 'hexane',\
        'decane', 'pentadecane', 'icosane'])
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        constants.Vcs = vc.tolist()
        constants.omegas = w.tolist()
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        kijs = Bin.tolist()

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)
        zs = [0.5, 0.03, 0.07, 0.2, 0.15, 0.05] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()
        PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo: ' +str(dt)+' s')

        t0_thermo = time.time()
        PT_vec = flasher.grid_flash(Ts=T, Ps=P, zs=z.T)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo vectorized: ' +str(dt)+' s')
        import pdb; pdb.set_trace()

    @unittest.skip("?")
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

    @unittest.skip("?")
    def test_caso2(self):
        R = 10.73159
        Bin = np.array([[0, 0.033],[0.033, 0]])
        Tc = np.array([369.8, 425.2])*9/5
        Pc = 14.7*np.array([41.9, 37.5])
        w = np.array([0.152, 0.193])
        Mw = np.array([44.097, 58.124])
        z = np.array([0.5, 0.5]) [:,np.newaxis]*np.ones([2,1])
        T = np.array([9/5*396])
        P = np.array([14.7*(3.86e6/101325)]) #* np.ones(10)

        print('\ncaso2:')
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P)
        obj.run(z,Mw)
        import pdb; pdb.set_trace()
        print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L,'V: ',obj.V)
        print('fl: ',obj.fl,'fv: ',obj.fv)

    @unittest.skip("ok")
    def test_caso3(self):
        R = 10.73159
        z = np.array([0.5,0.5])[:,np.newaxis]#*np.ones([2,10]) #exemplo aleatório
        P = np.array([(100*1E5)])/101325*14.7 #* np.ones(10)# pressão de 100bar e ele converte para atm e depois psi
        T = np.array([350])*9/5 #- T em K to R
        Tc = np.array([190.6, 460.4])*9/5;
        Pc =np.array([45.4,33.4])*14.7; # 14.7*class
        w = np.array([0.008,0.227])
        Bin = np.array([[0,0.0236],[0.0236,0]])
        Mw = np.array([44.097,58.124])


        print('\ncaso3:')
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P)
        obj.run(z,Mw)

        print('x: ',obj.x,'y: ',obj.y)
        print('K: ',obj.K)
        print('L: ',obj.L,'V: ',obj.V)
        print('fl: ',obj.fl,'fv: ',obj.fv)

    @unittest.skip("?")
    def test_bookDandekar(self):
         R = 10.73159
         z = np.array([0.8232,0.0871,0.0505,0.0198,0.0194])[:,np.newaxis]
         Tc = np.array([343,666,845,1112,1291])#Rankine
         Pc = np.array([667.2,615.8,489.4,305.7,205.7]) #psia
         Bin = np.array([[0,0.009,0.021,0.052,0.080],\
                         [0.009,0,0.003,0.019,0.039],\
                         [0.021,0.003,0,0.008,0.022],\
                         [0.052,0.019,0.008,0,0.004],\
                         [0.08,0.039,0.022,0.004,0]])
         w = np.array([0.008,0.152,0.251,0.490,0.742])
         Mw = np.array([16.043,44.097,72.15,142.29,226.41])
         T = np.array([100 + 459.67])#Fahrenheit
         P = np.array([1500]) #psia

         print('\ncasoN:')

         obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P)
         obj.run(z,Mw)
         '''
         obj.x = x
         obj.y = x
         obj.L = L #fração molar de liquido
         obj.V = V #fração molar de vapor
         obj.fl = fl #fugacidade dos componentes na fase liquida
         obj.fv = fv #fugacidade dos componentes na fase vapor
         '''
         print('x: ',obj.x,'y: ',obj.y)
         print('K: ',obj.K)
         print('L: ',obj.L,'V: ',obj.V)
         print('fl: ',obj.fl,'fv: ',obj.fv)
         import pdb; pdb.set_trace()

    @unittest.skip('ok')
    def test_case2_Moshiri(self):
        R = 8.3144598
        z = np.array([[0.9, 0.1, 0.0]]).T#[:,np.newaxis]#*np.ones([10]) #exemplo aleatório
        P = np.array([6.9e6])
        T = np.array([311])
        Tc = np.array([190.56, 305.32, 369.83])
        Pc = np.array([4599000, 4873000, 4248000])
        Mw = np.array([16.043e-3, 30.070e-3, 44.096e-3])
        w = np.array([0.0115, 0.0995, 0.1523])
        Bin = np.array([[0.,0.,.0], [0.,0.,0.], [0.,0.,.0]])
        print('\ncaso2 MM:')
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P)
        obj.run(z,Mw)
        import pdb; pdb.set_trace()

    @unittest.skip('ok')
    def test1_Moshiri(self):
        R = 8.3144598
        n = 10000
        z = np.array([[0.61048048, 0., 0.14185506, 0.13189786, 0.11576659]]).T# * np.ones((1,n))#[n_bXn_k]
        #P = np.linspace(100e6, 10000e6, n)
        P = np.array([20000e6]) #* np.ones(5)
        T = np.array([20+273])
        Tc = np.array([304.13, 190.56, 425.2, 618.1, 723])
        Pc = np.array([7380000, 4599000, 38e5, 21e5, 1.41e6])
        vc = np.array([0.0939e-3, 9.85988846e-05, 25.5e-5, 0.624e-3, 0.943e-3])
        Mw = np.array([44.01e-3,16.043e-3, 58.12e-3, 142.28e-3, 226.44e-3])
        w = np.array([0.225, 0.0115, 0.199, 0.484, 0.721])
        Bin = np.array([[0,0.,0.,.0, 0.],[0.,0.,0.,0.,0.],[0.,0.,.0,0.,0.], [0.,0.,0.,0.,0.], [0.,0.,0.,0.,0.]])
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P)
        obj.run(z,Mw)

        import pdb; pdb.set_trace()
        #[34121.18701394, 33506.58479809, 13289.01156059, 5171.6375957 , 2987.51943215])
        #[34433.22826536, 33892.84905336, 13342.42587675,  5180.17927084, 2990.4853912 ]
        #[34908.17670552, 34466.83197379, 13422.02659636,  5192.77066779, 2994.82879636]
        #[36032.03403002, 35757.3989117 , 13601.37719889,  5220.48095586, 3004.26222542]
    @unittest.skip('ok')
    def test_case1_BRB(self):
        R = 8.3144598
        z = np.array([0.95, 0.05])[:,np.newaxis]#*np.ones() #exemplo aleatório
        #P = np.array([4059356915.252799])
        P = np.array([25.65e6])
        #P = np.linspace(7.4e6, 30e6, 10)
        #T = np.array([25+273.15])
        T = np.array([299.82])
        Tc = np.array([304.21, 190.6])
        Pc = np.array([7.39e6, 4.6e6])
        Mw = np.array([44.01e-3, 16.04e-3])
        w = np.array([0.225, 0.022])
        Bin = np.array([[0., 0.12], [0.12, 0.]])
        vc = np.array([9.4e-5, 9.99e-05])

        print('\ncaso1 BRB:')
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P)
        obj.run(z,Mw)

        ''' Using Thermo '''
        constants, properties = ChemicalConstantsPackage.from_IDs(['CO2', 'methane'])
        constants.MWs = Mw.tolist()
        constants.Tcs = Tc.tolist()
        constants.Pcs = Pc.tolist()
        constants.Vcs = vc.tolist()
        constants.omegas = w.tolist()
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        kijs = Bin.tolist()

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)
        zs = [0.95, 0.05] #z.T.tolist()#[0.965, 0.018, 0.017]
        t0_thermo = time.time()

        PT = flasher.flash(T=T[0], P=P[0], zs=zs)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo: ' +str(dt)+' s')

        t0_thermo = time.time()
        PT_vec = flasher.grid_flash(Ts=T, Ps=P, zs=z.T)
        t1_thermo = time.time()
        dt = t1_thermo-t0_thermo
        print('Time for thermo vectorized: ' +str(dt)+' s')
        import pdb; pdb.set_trace()

        import pdb; pdb.set_trace()

    @unittest.skip("ok")
    def test_moyner(self):
        # CO2 C1 C10
        R = 8.3144598
        z = np.array([1, 0., 0.])[:,np.newaxis]
        Tc = np.array([304.21, 190.58, 617.65])
        Pc = np.array([7382235.525, 4603093.425, 2107154.7])
        vc = np.array([10.55286406e-5, 10.60259781e-05, 62e-5])
        Mw = np.array([44.429e-3, 16.043e-3, 142.9e-3])
        w = np.array([0.225, 0.008, 0.49])
        Bin = np.array([[0., 0., 0.], [0., 0., 0.0], [0., 0.0, .0]])
        P = np.array([125e5])
        T = np.array([423])
        print('\ncaso1:')
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P)
        obj.run(z,Mw)
        import pdb; pdb.set_trace()

    @unittest.skip("ok")
    def test_MY10(self):
        R = 8.3144598
        # C1, C2, C3, n-C4, n-C5, n-C6, n-C7, n-C8, n-C10, n-C14
        Nt = 100
        z = np.array([0.35, 0.03, 0.04, 0.06, 0.04, 0.03, 0.05, 0.05, 0.3, 0.05])[:,np.newaxis]
        Nk = z * Nt
        Tc = np.array([190.6, 305.4, 369.8, 425.2, 469.6, 507.5, 540.3, 568.8, 617.9, 691.9]) # Kelvin
        Pc = np.array([45.4, 48.2, 41.9, 37.5, 33.3, 30.1, 27.4, 24.9, 21.0, 15.2])*100000 # pascal
        Mw = np.array([16.04, 30.07, 44.1, 58.12, 72.15, 86.178, 100.205, 114.232, 142.29, 198.39])*1e-3
        w = np.array([0.008, 0.098, 0.152, 0.193, 0.251, 0.305, 0.305, 0.396, 0.484, 0.747])

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

        P = np.array([75.4])*100000
        T = np.array([566.6])

        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P)
        obj.run(z,Mw)
        x_or = obj.x
        y_or = obj.y
        L_or = obj.L
        V_or = obj.V
        Csi_l_or = obj.ksi_L
        Csi_v_or = obj.ksi_V
        K_or = obj.K
        import pdb; pdb.set_trace()
        delta = 0.001
        Nk_plus = np.copy(Nk)
        Nk_minus = np.copy(Nk)
        for i in range(len(z)):
            Nk_plus[i,0] = Nk[i,0] + delta/2
            Nk_minus[i,0] = Nk[i,0] - delta/2
            z_plus = Nk_plus/sum(Nk_plus)
            z_minus = Nk_minus/sum(Nk_minus)
            obj.run(z_plus,Mw)
            x_plus = obj.x
            y_plus = obj.y
            L_plus = obj.L
            V_plus = obj.V
            Csi_l_plus = obj.ksi_L
            Csi_v_plus = obj.ksi_V
            K_plus = obj.K
            import pdb; pdb.set_trace()


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

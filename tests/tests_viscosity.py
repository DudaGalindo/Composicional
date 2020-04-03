import numpy as np
from viscosity import LorenzBrayClark
from stability_check import StabilityCheck
import unittest

class Tests_Viscosity(unittest.TestCase):
    def test1_methane_propane(self):
        R = 8.3144598
        z = np.array([0.2, 0.8])
        T = 25 + 273
        P = 69.1 * 101325
        Mw = np.array([16.0425, 44.1]) *10**(-3)
        Pc = np.array([46, 42.5]) * 1e5
        Tc = np.array([190.6, 370])
        w = np.array([0.011, 0.153])
        vc = np.array([99.6, 203]) * 10**(-6)
        Bin = np.zeros([len(w),len(w)])
        C7 = np.array([0, 0])

        print('\nmet_prop:\n')
        L = 1; V = 0
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P,C7)
        obj.run(z,Mw)
        #obj.run_other_properties(z, L, V, Mw) #caso que é líquido e só. Não teria check de estabilidade

        component_molar_fractions = np.zeros([len(w),2,1])
        component_molar_fractions[:,0,0] = obj.x
        component_molar_fractions[:,1,0] = obj.y


        #obj.rho_L = 420.5 # tem que acrescentar vshift
        #obj.x = z
        #obj.ksi_L = obj.rho_L / sum(obj.x * Mw)

        phase_molar_densities = np.zeros([1,2,1])
        phase_molar_densities[0,0,:] = obj.ksi_L
        phase_molar_densities[0,1,:] = obj.ksi_V


        phase_viscosities = LorenzBrayClark(Mw, Pc, vc, Tc, P, T, phase_molar_densities, component_molar_fractions)
        mis = phase_viscosities(obj)
        mi = obj.L * mis[0,0,0] + obj.V * mis[0,1,0]
        import pdb; pdb.set_trace()

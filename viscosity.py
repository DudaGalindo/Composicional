import numpy as np


class LorenzBrayClark:

    def __init__(self, Mw, Pc, vc, Tc, P, T, phase_molar_densities, component_molar_fractions):
        self.n_volumes = 1
        self.phase_molar_densities = phase_molar_densities * 10**(-6) # mole/m³ to mole/cm³
        self.Mw =  Mw * 10**(3)  # Kg/mole to grams/mole
        self.Pc = Pc / 101325 # Pa to atm
        self.vc = vc * 10 ** 6 # m³/mole to cm³/mole
        self.Tc = Tc
        self.P = P / 101325
        self.T = T
        self.Nc = len(Mw)
        self.component_molar_fractions = component_molar_fractions

    def component_viscosity(self):
        mi_components = np.zeros(self.Nc)
        Trs = self.T / self.Tc
        ind_Tr_lower = np.argwhere(Trs <= 1.5)
        ind_Tr_higher = np.argwhere(Trs > 1.5)

        # in this formula, the component molecular weight is in grams/mole. and
        # the pressure is in atm.
        # self.Mw and self.Pc are in kg/mole and Pa, respectively.
        pure_components_viscosity_parameters = self.Tc ** (1/6) \
                 / ((self.Mw) ** (1/2) * (self.Pc) ** (2/3))

        mi_components[ind_Tr_lower] = 3.4e-4 * Trs[ind_Tr_lower] ** 0.94 / \
                        pure_components_viscosity_parameters[ind_Tr_lower]
        mi_components[ind_Tr_higher] = (1.778e-4 * (4.58 * Trs[ind_Tr_higher] -
        1.67) ** (5/8)) / pure_components_viscosity_parameters[ind_Tr_higher]

        self.mi_components = mi_components[:, np.newaxis, np.newaxis]
        #self.mi_components = np.ones([self.Nc, 1, self.n_volumes]) * mi_components

    def phase_viscosity_atm(self, prop):

        self.mi_atm = np.sum(self.component_molar_fractions * self.mi_components *
                self.Mw[:, np.newaxis, np.newaxis] ** (1/2), axis = 0) \
                /np.sum(self.component_molar_fractions *
                self.Mw[:, np.newaxis, np.newaxis] ** (1/2) , axis = 0)

        self.mi_atm = self.mi_atm[np.newaxis,:,:]

    def phase_viscosity(self):
        #include in the entry parameters the vc: component critical molar volume
        # mi_phase = np.zeros([1,2,self.n_volumes])
        a = np.array([0.1023, 0.023364, 0.058533, -0.040758, 0.0093324])
        #self.phase_molar_densities = np.zeros([1, 2, self.n_volumes])
        #self.phase_molar_densities[0,0,:] = self.phase_molar_densities[:,0,:]
        #self.phase_molar_densities[0,1,:] = self.phase_molar_densities[:,1,:]

        phase_reduced_molar_density = self.phase_molar_densities[:,0,:] * \
            np.sum(self.component_molar_fractions * self.vc[:,np.newaxis,np.newaxis], axis = 0)

        # ind_lower = np.argwhere(phase_reduced_molar_density <= 0.18)
        # ind_higher = np.argwhere(phase_reduced_molar_density > 0.18)

        Xs = a[0]+ a[1] * phase_reduced_molar_density + \
            a[2] * phase_reduced_molar_density ** 2 + \
            a[3] * phase_reduced_molar_density ** 3 + \
            a[4] * phase_reduced_molar_density ** 4

        # parametro de viscosidade para mistura

        neta = (np.sum(self.component_molar_fractions *
                self.Tc[:,np.newaxis,np.newaxis] , axis = 0)) ** (1/6) \
                / (np.sum(self.component_molar_fractions *
                self.Mw[:,np.newaxis,np.newaxis], axis = 0) ** (1/2)
                * np.sum(self.component_molar_fractions *
                self.Pc[:,np.newaxis,np.newaxis], axis = 0) ** (2/3))
        # mi_phase[ind_lower] = self.mi_mix[ind_lower] + 2.05e-4 * \
        #                             self.zetas_r[ind_lower] / neta[ind_lower]
        # mi_phase[ind_higher] = self.mi_mix[ind_higher] + (Xs[ind_higher] ** 4 - 1) \
        #                         /(1e4*neta[ind_higher])

        mi_phase = (self.mi_atm + (Xs ** 4 - 1e-4) / neta) * 1e-3 #return in Pa.s
        import pdb; pdb.set_trace()
        return mi_phase

    def __call__(self, prop):
        self.component_viscosity()
        self.phase_viscosity_atm(prop)
        mi_phase = self.phase_viscosity()
        return mi_phase

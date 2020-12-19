"""Check stability of a thermodynamic equilibrium."""
import numpy as np
import math
import thermo
from scipy.misc import derivative
from equation_of_state import PengRobinson
import matplotlib.pyplot as plt
## Encontrar os pontos estacionarios. Estes correspondem aos pontos nos quais a derivada de g com respeito a Y é 0
## Todas as equações foram confirmadas pelo livro do Dandekar e a biblioteca do thermo

class StabilityCheck:
    """Check for stability of a thermodynamic equilibrium and returns the
    equilibrium phase compositions (perform the flash calculation)."""

    def __init__(self, w, Bin, R, Tc, Pc, T, P, Pb_guess):
        self.w = w
        self.Bin = Bin
        self.R = R
        self.Tc = Tc
        self.Pc = Pc
        self.T = T
        self.P = P
        self.Nc = len(w)
        #self.Pv = Pv
        self.ph_L = np.ones(len(self.P), dtype = bool)
        self.ph_V = np.zeros(len(self.P), dtype = bool)
        self.Pb_guess = Pb_guess


        #StabilityCheck.TPD(self)

    def run(self, z, Mw):
        PR = PengRobinson(self)
        self.vapor_pressure_pure_substancies(PR)
        self.equilibrium_ratio_Wilson()
        self.L = np.empty(len(self.P))
        self.V = np.empty(len(self.P))
        self.x = np.empty(z.shape)
        self.y = np.empty(z.shape)

        ponteiro_flash = np.zeros(len(self.P), dtype = bool)
        dir_flash = np.argwhere(z <= 0)
        ponteiro_flash[dir_flash[:,1]] = True


        if any(~ponteiro_flash):
            sp1,sp2 = self.Stability(PR, z, np.copy(~ponteiro_flash))
            ponteiro_aux = ponteiro_flash[~ponteiro_flash]
            ponteiro_aux[(np.round(sp1,14) > 1) + (np.round(sp2,14) > 1)] = True #os que devem passar para o calculo de flash
            ponteiro_flash[~ponteiro_flash] = ponteiro_aux

        self.molar_properties(PR, z, np.copy(ponteiro_flash))

        ponteiro_flash[self.L<0] = False
        ponteiro_flash[self.L>1] = False
        import pdb; pdb.set_trace()
        self.x[:,~ponteiro_flash] = z[:,~ponteiro_flash]
        self.y[:,~ponteiro_flash] = z[:,~ponteiro_flash]
        self.bubble_point_pressure(PR, z, Mw, np.copy(~ponteiro_flash))
        #self.x = z
        #self.y = z
        #self.L = 1
        #self.V = 0
        self.get_other_properties(PR, Mw)
        import pdb; pdb.set_trace()
        #PR.get_all_derivatives(self, self.x, self.y, self.L, self.V, self.P)

    def get_other_properties(self, PR, Mw):
        self.Mw_L, self.ksi_L, self.rho_L = self.other_properties(PR, self.x, Mw, self.ph_L)
        self.Mw_V, self.ksi_V, self.rho_V = self.other_properties(PR, self.y, Mw, self.ph_V)

    def equilibrium_ratio_Wilson(self):

        self.K = np.exp(5.37 * (1 + self.w) * (1 - 1 / (self.T / self.Tc)), dtype=np.double)[:,np.newaxis] / \
                (self.P / self.Pc[:,np.newaxis])


    """-------------Below starts stability test calculation -----------------"""

    def Stability(self, PR, z, ponteiro_stab_check):
        ''' In the lnphi function: 0 stands for vapor phase and 1 for liquid '''
    #****************************INITIAL GUESS******************************#
    ## Both approaches bellow should be used in case the phase is in the critical region

    #*****************************Test one**********************************#
        #Used alone when the phase investigated (z) is clearly vapor like (ph = 0)
        Y = np.empty(z.shape)
        lnphiz = np.empty(z.shape)
        ponteiro = np.copy(ponteiro_stab_check)
        Y[:,ponteiro] = z[:,ponteiro] / self.K[:,ponteiro]
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphiz[:,ponteiro] = PR.lnphi(self, z[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = PR.lnphi(self, y[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(z[:,ponteiro]) + lnphiz[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        stationary_point1 = np.sum(Y[:,ponteiro_stab_check], axis = 0)


        # print(sum(Y))
        # if sum(Y) <= 1: print('estavel')
        # else: print('instavel') #estavel2 = 0

        ''' If this first test returns stable: There might be a chance that the
        Gibbs free energy is at its global minimum. However, as will be
        explaned later, this alone doesn't guarantee the phase stability.
         If it returns unstable: There is another composition that makes the
        Gibbs free energy at its global minimum indicated by sum(Y)>1'''
    #*****************************Test two**********************************#
        #Used alone when the phase investigated (z) is clearly liquid like (ph == 1)
        ponteiro = np.copy(ponteiro_stab_check)

        Y[:,ponteiro] = self.K[:,ponteiro] * z[:,ponteiro]
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphiz[:,ponteiro] = PR.lnphi(self, z[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = PR.lnphi(self, y[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(z[:,ponteiro]) + lnphiz[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        stationary_point2 = np.sum(Y[:,ponteiro_stab_check], axis = 0)

        L = self.L[ponteiro_stab_check]
        L[(np.round(stationary_point1,4)==1) * (np.round(stationary_point2,4)<=1)] = 1.
        L[(np.round(stationary_point1,4)<1) * (np.round(stationary_point2,4)==1)] = 0.
        self.L[ponteiro_stab_check] = L
        self.V[ponteiro_stab_check] = 1 - self.L[ponteiro_stab_check]
        self.x[:,ponteiro_stab_check] = z[:,ponteiro_stab_check]
        self.y[:,ponteiro_stab_check] = z[:,ponteiro_stab_check]
        return stationary_point1, stationary_point2

        """
        The same thing happens here. The difference is that, the original
        phase is gas, and then the "new" phase is supposed to be liquid.
        In cases that the compressibility equation returns only one root,
        both tests work like two different initial guess for the same problem,
        being more likely to find a stationary point that makes the phase
        unstable.
        """

        ''' If one of these approaches returns unstable the system is unstable.
        The stability of the phase is something a little bit more complex
        to guarantee. '''
        #if stationary_point1 == 1: return stationary_point1,stationary_point2, Pb1
        #if stationary_point2 == 1: return stationary_point1,stationary_point2, Pb2
        return stationary_point1,stationary_point2

    """-------------Below starts biphasic flash calculations-----------------"""
    def molar_properties(self, PR, z, ponteiro):
        if self.Nc <= 7: ponteiro = self.molar_properties_Whitson(PR, z, ponteiro)
        else: ponteiro = self.molar_properties_Yinghui(PR, z, ponteiro)
        return ponteiro

    def deltaG_molar_vectorized(self, PR, l, P, ph):
        lnphi = np.empty([2, len(self.w), len(ph)])

        lnphi[0,:] = PR.lnphi(self, l, P, 1 - ph)
        lnphi[1,:] = PR.lnphi(self, l, P, ph)
        deltaG_molar = np.sum(l * (lnphi[1 - ph,:, np.arange(len(ph))] - lnphi[1*ph,:, np.arange(len(ph))]).T, axis = 0)
        ph[deltaG_molar<0] = 1 - ph[deltaG_molar<0]
        """
        This part above means that:
        if deltaGx_molar>=0: lnphil = lnphix[1]
                            (lnphi for x in the liquid phase)
        if deltaGx_molar<0: lnphil = lnphix[0]
                            (lnphi for x in the vapour phase)
                            because it minimizes the Gibbs molar energy
        The same analysis is made for the y phase composition.
        """
        return ph


    def lnphi_based_on_deltaG(self, PR, l, P, ph):
        ph = self.deltaG_molar_vectorized(PR, l, P, ph)
        return PR.lnphi(self, l, P, ph)

    def solve_objective_function_Yinghui(self, z1, zi, K1, KNc, Ki, K, x):
        x1_min = z1 * (1 - KNc) / (K1 - KNc)
        x1_max = (1 - KNc) / (K1 - KNc)

        vols_zi_neg = np.zeros(len(K1), dtype = bool)
        vols_zi_neg[np.sum(zi < 0, axis = 0, dtype=bool)] = True
        KNc_z_neg = KNc[vols_zi_neg]
        K1_z_neg = K1[vols_zi_neg]
        z1_z_neg = z1[vols_zi_neg]
        zi_z_neg = zi[:,vols_zi_neg]
        Ki_z_neg = Ki[:,vols_zi_neg]
        vols_zi_neg_num = np.sum(vols_zi_neg*1) + 1 - np.sign(np.sum(vols_zi_neg*1))
        Ki_z_neg_K_big1 = Ki_z_neg[Ki_z_neg > 1].reshape(int(len(Ki_z_neg[Ki_z_neg>1])/vols_zi_neg_num), vols_zi_neg_num)
        theta = np.ones(zi[:,vols_zi_neg].shape)

        theta[Ki_z_neg > 1] = ((1 - KNc_z_neg[np.newaxis,:]) / (Ki_z_neg_K_big1 - KNc_z_neg[np.newaxis,:])).ravel()
        aux_eq = (Ki_z_neg - 1) * z1_z_neg[np.newaxis,:] / (zi_z_neg * (K1_z_neg[np.newaxis,:] - 1) /
                theta - (K1_z_neg[np.newaxis,:] - Ki_z_neg))

        #aux_eq = (K - 1) * z1 / (z * (K1 - 1) / theta - (K1 - K))
        cond = (Ki_z_neg[zi_z_neg != 0] - 1) * z1_z_neg[np.newaxis,:] / zi_z_neg[zi_z_neg != 0]
        cond_aux = np.ones(cond.shape[1], dtype = bool)
        cond_aux[np.sum(cond <= 0, axis = 0, dtype=bool)] = False
        aux_eq_cond = aux_eq[cond_aux]
        vols_aux = len(cond_aux==True)

        vols_aux = np.sum(cond_aux*1) + 1 - np.sign(np.sum(cond_aux*1))
        aux_eq_cond = aux_eq_cond[aux_eq_cond >= 0].reshape(int(len(aux_eq_cond[aux_eq_cond >= 0])/vols_aux), vols_aux)
        x1_max_aux = np.copy(x1_max[vols_zi_neg])
        if any(cond_aux): x1_max_aux[cond_aux] = np.min(aux_eq_cond, axis = 0)
        x1_max[vols_zi_neg] = x1_max_aux

        x1_min_aux = np.copy(x1_min[vols_zi_neg])
        if any(~cond_aux): x1_min_aux[~cond_aux] = np.max(aux_eq[~cond_aux], axis = 0)
        x1_min_aux[x1_min_aux < 0] = 0
        x1_min[vols_zi_neg] = x1_min_aux

        if any(x1_min > x1_max): raise ValueError('There is no physical root')

        x1 = (x1_min + x1_max) / 2

        ponteiro = np.ones(len(x1), dtype = bool)

        while any(ponteiro):
            f = 1 + ((K1[ponteiro] - KNc[ponteiro]) / (KNc[ponteiro] - 1)) * x1[ponteiro] + np.sum(((Ki[:,ponteiro] - KNc[ponteiro][np.newaxis,:]) /
                (KNc[ponteiro][np.newaxis,:] - 1)) * zi[:,ponteiro] * (K1[ponteiro][np.newaxis,:] - 1) * x1[ponteiro][np.newaxis,:]
                / ((Ki[:,ponteiro] - 1) * z1[ponteiro][np.newaxis,:] + (K1[ponteiro][np.newaxis,:] - Ki[:,ponteiro]) *
                x1[ponteiro][np.newaxis,:]), axis = 0)
            df = ((K1[ponteiro] - KNc[ponteiro]) / (KNc[ponteiro] - 1)) + np.sum(((Ki[:,ponteiro] - KNc[ponteiro][np.newaxis,:]) /
                (KNc[ponteiro][np.newaxis,:] - 1)) * zi[:,ponteiro] * z1[ponteiro][np.newaxis,:] * (K1[ponteiro][np.newaxis,:] - 1) *
                (Ki[:,ponteiro] - 1) / ((Ki[:,ponteiro] - 1) * z1[ponteiro][np.newaxis,:] + (K1[ponteiro][np.newaxis,:] - Ki[:,ponteiro]) *
                x1[ponteiro][np.newaxis,:]) ** 2, axis = 0)
            x1[ponteiro] = x1[ponteiro] - f/df #Newton-Raphson iterative method
            x1_aux = x1[ponteiro]
            x1_aux[x1_aux > x1_max] = (x1_min[x1_aux > x1_max] + x1_max[x1_aux > x1_max])/2
            x1_aux[x1_aux < x1_min] = (x1_min[x1_aux < x1_min] + x1_max[x1_aux < x1_min])/2
            x1[ponteiro] = x1_aux
            ponteiro_aux = ponteiro[ponteiro] #o que muda de tamanho
            ponteiro_aux[f < 1e-10] = False
            ponteiro[ponteiro] = ponteiro_aux
            x1_max = x1_max[ponteiro_aux]
            x1_min = x1_min[ponteiro_aux]
            x1_max[f[ponteiro_aux] * df[ponteiro_aux] > 0] = x1[ponteiro][f[ponteiro_aux] * df[ponteiro_aux] > 0]
            x1_min[f[ponteiro_aux] * df[ponteiro_aux] < 0] = x1[ponteiro][f[ponteiro_aux] * df[ponteiro_aux] < 0]


        xi = (K1[np.newaxis,:] - 1) * zi * x1[np.newaxis,:] / ((Ki - 1) * z1[np.newaxis,:] +
            (K1[np.newaxis,:] - Ki) * x1[np.newaxis,:])

        x_not_z1_zero = np.copy(x)
        x_not_z1_zero[K == K1[np.newaxis,:]] = x1
        x_not_z1_zero[K == KNc[np.newaxis,:]] = 1 - np.sum(xi, axis = 0) - x1
        aux_xi = np.ones(x_not_z1_zero.shape,dtype=bool)
        aux_xi[K == K1[np.newaxis,:]] = False
        aux_xi[K == KNc[np.newaxis,:]] = False
        x_not_z1_zero[aux_xi] = xi.ravel()
        return x_not_z1_zero

    def Yinghui_method(self, z, ponteiro):

        """ Shaping K to Nc-2 components by removing K1 and KNc and z to Nc-2
        components by removing z1 and zNc """
        K = self.K[:,ponteiro]
        x = self.x[:,ponteiro]
        K1 = np.max(K, axis = 0); KNc = np.min(K, axis=0)
        z1 = z[K == K1[np.newaxis,:]]

        aux = np.ones(K.shape, dtype = bool)
        aux[K == K1[np.newaxis,:]] = False
        aux[K == KNc[np.newaxis,:]] = False
        Ki = K[aux]
        zi = z[aux]


        vols_ponteiro = np.sum(ponteiro*1) + 1 - np.sign(np.sum(ponteiro*1))
        Ki = Ki.reshape(int(len(Ki)/vols_ponteiro), vols_ponteiro)
        zi = zi.reshape(int(len(zi)/vols_ponteiro), vols_ponteiro)

        #starting x

        x[:,~(z1 == 0)] = self.solve_objective_function_Yinghui(z1[~(z1 == 0)], zi[:,~(z1 == 0)], K1[~(z1 == 0)], KNc[~(z1 == 0)], Ki[:,~(z1 == 0)],
                                            K[:,~(z1 == 0)], x[:,~(z1 == 0)])


        #self.solve_objective_function_Yinghui_explicitly()
        z_z1_zero = z[:,z1 == 0]
        K_z1_zero = K[:,z1 == 0]
        K_KNc_z1_zero = K_z1_zero[K_z1_zero == KNc[z1 == 0][np.newaxis,:]]

        aux_xNc = np.zeros(K_z1_zero.shape, dtype = bool); aux_x1 = np.copy(aux_xNc)
        aux_xNc[K_z1_zero == KNc[z1 == 0][np.newaxis,:]] = True
        aux_x1[K_z1_zero == K1[z1 == 0][np.newaxis,:]] = True
        aux_xi = ~(aux_xNc + aux_x1)
        xi_z1_zero = ((K1[z1 == 0][np.newaxis,:] - 1) * zi[:,z1 == 0] / (K1[z1 == 0][np.newaxis,:] - Ki[:,z1 == 0]))
        x_z1_zero = np.zeros(x[:,z1 == 0].shape)
        x_z1_zero[aux_xNc] = (K1[z1 == 0] - 1) * z_z1_zero[aux_xNc] / (K1[z1 == 0] - K_z1_zero[aux_xNc])
        x_z1_zero[aux_xi] = xi_z1_zero.ravel()
        x_z1_zero[aux_x1] = 1 - np.sum(x_z1_zero, axis = 0)
        x[:,z1 == 0] = x_z1_zero
        self.x[:,ponteiro] = x
        self.y[:,ponteiro] = self.K[:,ponteiro] * self.x[:,ponteiro]

    def molar_properties_Yinghui(self, PR, z, ponteiro):
        #razao = fl/fv -> an arbitrary vector to enter in the iterative mode

        razao = np.ones(z.shape)/2
        ponteiro_save = np.copy(ponteiro)

        while any(ponteiro):
            self.Yinghui_method(z[:,ponteiro], ponteiro)
            lnphil = self.lnphi_based_on_deltaG(PR, self.x[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            lnphiv = self.lnphi_based_on_deltaG(PR, self.y[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            self.fl = np.exp(lnphil) * (self.x[:,ponteiro] * self.P[ponteiro][np.newaxis,:])
            self.fv = np.exp(lnphiv) * (self.y[:,ponteiro] * self.P[ponteiro][np.newaxis,:])
            razao[:,ponteiro] = np.divide(self.fl, self.fv, out = razao[:,ponteiro] / razao[:,ponteiro] * (1 + 1e-10),
                              where = self.fv != 0)
            self.K[:,ponteiro] = razao[:,ponteiro] * self.K[:,ponteiro]
            stop_criteria = np.max(abs(razao[:,ponteiro] - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        V = (z[:,ponteiro_save][self.x[:,ponteiro_save] != 0] - self.x[:,ponteiro_save][self.x[:,ponteiro_save] != 0]) / \
                          (self.y[:,ponteiro_save][self.x[:,ponteiro_save] != 0] - self.x[:,ponteiro_save][self.x[:,ponteiro_save] != 0])
        #vols_V = np.sum(self.x[:,ponteiro_save] == 0, dtype=bool, axis = 0)
        vv = np.argwhere((self.x[:,ponteiro_save]!=0) == True)
        vols_V, ind = np.unique(vv[:,1],return_index = True)
        self.V[ponteiro_save] = V[ind]
        self.L[ponteiro_save] = 1. - self.V[ponteiro_save]
        return ponteiro_save

    def solve_Whitson_for_L(self, z, ponteiro):
        K = self.K[:,ponteiro]
        Lmax = np.max(K, axis = 0)/(np.max(K, axis = 0) - 1)
        Lmin = np.min(K, axis = 0)/(np.min(K, axis = 0) - 1)
        L = (Lmin + Lmax) / 2

        ponteiro_save = np.copy(ponteiro)
        while any(ponteiro):
        #while abs(L / Lold - 1) > 1e-9:
             Lold = L[ponteiro]
             f = np.sum((1 - K[:,ponteiro]) * z[:,ponteiro] / (L[ponteiro][np.newaxis,:]
                        + (1 - L[ponteiro][np.newaxis,:]) * K[:,ponteiro]), axis = 0)
             df = - np.sum((1 - K[:,ponteiro]) ** 2 * z[:,ponteiro] / (L[ponteiro][np.newaxis,:] +
                    (1 - L[ponteiro][np.newaxis,:]) * K[:,ponteiro]) ** 2, axis = 0)
             L[ponteiro] = L[ponteiro] - f / df #Newton-Raphson iterative method
             L_aux = L[ponteiro]
             L_aux[L_aux > Lmax[ponteiro]] = (Lmax[ponteiro] + L[ponteiro]) / 2
             L_aux[L_aux < Lmin[ponteiro]] = (Lmin[ponteiro] + L[ponteiro]) / 2
             L[ponteiro] = L_aux
             stop_criteria = abs(L[ponteiro] / Lold - 1)
             ponteiro_aux = ponteiro[ponteiro]
             ponteiro_aux[stop_criteria < 1e-9] = False
             ponteiro[ponteiro] = ponteiro_aux
        self.L[ponteiro_save] = L[ponteiro_save]
        self.V[ponteiro_save] = (1 - self.L[ponteiro_save])
        self.x[:,ponteiro_save] = z[:,ponteiro_save] / (self.L[ponteiro_save][np.newaxis,:] +
                        (1 - self.L[ponteiro_save][np.newaxis,:]) * self.K[:,ponteiro_save])
        self.y[:,ponteiro_save] = self.K[:,ponteiro_save] * self.x[:,ponteiro_save]


    def solve_objective_function_Whitson_for_V(self, z, V, Vmax, Vmin, ponteiro):

        ponteiro_save = np.copy(ponteiro)
        while any(ponteiro):
            Vold = V[ponteiro]
            f = np.sum((self.K[:,ponteiro] - 1) * z[:,ponteiro] / (1 + V[ponteiro][np.newaxis,:] * (self.K[:,ponteiro] - 1)), axis = 0)
            df = - np.sum((self.K[:,ponteiro] - 1) ** 2 * z[:,ponteiro] / (1 + V[ponteiro][np.newaxis,:] * (self.K[:,ponteiro] - 1)) ** 2, axis = 0)
            self.V[ponteiro] = self.V[ponteiro] - f / df #Newton-Raphson iterative method
            V_aux = self.V[ponteiro]
            V_aux[V_aux > Vmax[ponteiro]] = 0.5 * (Vmax[ponteiro][V_aux > Vmax[ponteiro]] + Vold[V_aux > Vmax[ponteiro]]) #(Vmax + Vold)/2
            V_aux[V_aux < Vmin[ponteiro]] = 0.5 * (Vmin[ponteiro][V_aux < Vmin[ponteiro]] + Vold[V_aux < Vmin[ponteiro]])#(Vmax + Vold)/2
            self.V[ponteiro] = V_aux
            stop_criteria = abs(V[ponteiro] / Vold - 1)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux
            #import pdb; pdb.set_trace()
        #self.V[ponteiro_save] = V[ponteiro_save] + 1e-15 #manipulação
        self.x[:,ponteiro_save] = z[:,ponteiro_save] / (1 + self.V[ponteiro_save][np.newaxis,:] * (self.K[:,ponteiro_save] - 1))
        self.y[:,ponteiro_save] = self.K[:,ponteiro_save] * self.x[:,ponteiro_save]
        self.L[ponteiro_save] = 1. - self.V[ponteiro_save]


    def molar_properties_Whitson(self, PR, z, ponteiro):
        Lmax = np.max(self.K, axis = 0)/(np.max(self.K, axis = 0) - 1)
        Lmin = np.min(self.K, axis = 0)/(np.min(self.K, axis = 0) - 1)
        Vmax = 1. - Lmin
        Vmin = 1. - Lmax
        #Vmin = ((K1-KNc)*z[self.K==K1]-(1-KNc))/((1-KNc)*(K1-1))
        #proposed by Li et al for Whitson method
        self.V[ponteiro] = (Vmin[ponteiro] + Vmax[ponteiro]) * 0.5

        razao = np.ones(z.shape)/2
        ponteiro_save = np.copy(ponteiro)
        while any(ponteiro):
            self.solve_objective_function_Whitson_for_V(z, self.V, Vmax, Vmin, np.copy(ponteiro))
            lnphil = self.lnphi_based_on_deltaG(PR, self.x[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            lnphiv = self.lnphi_based_on_deltaG(PR, self.y[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            self.fv = np.exp(lnphiv) * (self.y[:,ponteiro] * self.P[ponteiro][np.newaxis,:])
            self.fl = np.exp(lnphil) * (self.x[:,ponteiro] * self.P[ponteiro][np.newaxis,:])
            razao[:,ponteiro] = np.divide(self.fl, self.fv, out = razao[:,ponteiro] / razao[:,ponteiro] * (1 + 1e-10),
                              where = self.fv != 0)
            self.K[:,ponteiro] = razao[:,ponteiro] * self.K[:,ponteiro]
            stop_criteria = np.max(abs(self.fv/self.fl - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux
            ponteiro[abs(self.L) > 2] = False

        return ponteiro_save

    def other_properties(self, PR, l, Mw, ph):
        #l - any phase molar composition

        A, B = PR.coefficients_cubic_EOS_all(self, l, self.P)
        ph = self.deltaG_molar_vectorized(PR, l, self.P, ph)
        Zfunc = np.vectorize(PengRobinson.Z)
        Z = Zfunc(B, A, ph)
        ksi_phase = self.P / (Z * self.R * self.T)
        Mw_phase = np.sum(l * Mw[:,np.newaxis], axis=0)
        rho_phase = ksi_phase * Mw_phase

        """
        Legenda:
            rho: phase mass density
            Mw: phase molecular weight
            eta: phase molar density
        """
        # se precisar retornar mais coisa, entra aqui
        return Mw_phase, ksi_phase, rho_phase

    def vapor_pressure_pure_substancies(self, PR):
        '''Lee-Kesler Correlation - only valid for T < Tc'''

        Tr = self.T/self.Tc
        #Tr[Tr>1] = 1

        A = 5.92714 - 6.09648 / Tr - 1.2886 * np.log(Tr) + 0.16934 * Tr**6
        B = 15.2518 - 15.6875 / Tr - 13.4721 * np.log(Tr) + 0.4357 * Tr**6
        self.Pv = self.Pc * np.exp(A + self.w * B)
        #self.Pv[self.T<self.Tc] = Pv[self.T<self.Tc]

        PR_kC7 = np.array([0.379642, 1.48503, 0.1644, 0.016667])
        PR_k = np.array([0.37464, 1.54226, 0.26992])


        '''k = (PR_kC7[0] + PR_kC7[1] * self.w - PR_kC7[2] * self.w ** 2 + \
            PR_kC7[3] * self.w ** 3) * (1*(self.w >= 0.49)) + (PR_k[0] + PR_k[1] * self.w - \
            PR_k[2] * self.w ** 2) * (1*(self.w < 0.49))
        alpha = (1 + k * (1 - (self.T / self.Tc) ** (1 / 2))) ** 2
        A_B = 0.45724/0.07780 * alpha / (self.T / self.Tc)
        razao = np.ones((self.Nc,1))
        xkj = np.zeros((self.Nc,1))

        for i in range(self.Nc):
            B[i] = 1.
            stop_criteria = 1
            xkj[i] = 1
            Pv[i] = 1.01618165e+09
            while stop_criteria > 1e-9:
                A, B = PR.coefficients_cubic_EOS_all(self, xkj, Pv[i])
                Zfunc = np.vectorize(PengRobinson.Z)
                Zl = Zfunc(B, A, np.array([1]))
                Zv = Zfunc(B, A, np.array([0]))
                import pdb; pdb.set_trace()
                lnphil = PR.lnphi_calculation(A, np.array([B[i]]), Zl)[i]
                lnphiv = PR.lnphi_calculation(A, np.array([B[i]]), Zv)[i]
                fil = np.exp(lnphil) * Pv[i]
                fiv = np.exp(lnphiv) * Pv[i]
                razao = np.divide(fil, fiv, out = razao / razao * (1 + 1e-10),
                                  where = fiv != 0)
                Pv = fil/fiv * Pv
                stop_criteria = np.max(abs(fiv/fil - 1), axis = 0)
            xkj[i] = 0'''

    def get_dlnphidP(self, PR, xij, P, ph):
        A, B = PR.coefficients_cubic_EOS_all(self, xij, P)
        Zfunc = np.vectorize(PengRobinson.Z)
        #ph = self.deltaG_molar_vectorized(PR, xij, P, ph)
        Z = Zfunc(B, A, ph)
        dAdP = PR.dA_dP(self)
        dBdP = PR.dB_dP(self)
        dZdP = PR.dZ_dP_parcial(dAdP, dBdP, Z, A, B)
        dlnphidP = PR.dlnphi_dP(self, dAdP, dBdP, dZdP, Z, A, B)
        return dlnphidP

    def bubble_point_pressure(self, PR, z, Pvi, ponteiro):
        ponteiro_save = np.copy(ponteiro)
        self.x[:,ponteiro] = z[:,ponteiro]#self.x[:,ponteiro]
        y = np.copy(self.x)
        # Depende muito de Pbguess
        Pb = self.Pb_guess*np.ones(len(self.P))

        K = np.exp(5.37 * (1 + self.w) * (1 - 1 / (self.T / self.Tc)), dtype=np.double)[:,np.newaxis] / \
            (Pb / self.Pc[:,np.newaxis])
        i = 0

        #Y[:,ponteiro] = K[:,ponteiro] * z[:,ponteiro]

        '''y[:,ponteiro] = K[:,ponteiro] * z[:,ponteiro]
        while any(ponteiro):
            lnphiv = self.lnphi_based_on_deltaG(PR, y[:,ponteiro], Pb[ponteiro], self.ph_V[ponteiro])
            lnphil = self.lnphi_based_on_deltaG(PR, self.x[:,ponteiro], Pb[ponteiro], self.ph_L[ponteiro])
            fil = np.exp(lnphil) * (self.x[:,ponteiro] * Pb[ponteiro][np.newaxis,:])
            fiv = np.exp(lnphiv) * (y[:,ponteiro] * Pb[ponteiro][np.newaxis,:])
            phiv = np.exp(lnphiv)
            phil = np.exp(lnphil)
            Pb[ponteiro] = np.sum(fil/phiv, axis = 0)
            y[:,ponteiro] = 1/Pb[ponteiro] * fil/phiv
            stop_criteria = np.sum((1 - fil/fiv)**2, axis=0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria <= 1e-10] = False
            ponteiro[ponteiro] = ponteiro_aux'''


        while any(ponteiro):

            y[:,ponteiro] = self.x[:,ponteiro] * K[:,ponteiro]
            Pb_old = np.copy(Pb[ponteiro])
            lnphiv = self.lnphi_based_on_deltaG(PR, y[:,ponteiro], Pb[ponteiro], self.ph_V[ponteiro])
            lnphil = self.lnphi_based_on_deltaG(PR, self.x[:,ponteiro], Pb[ponteiro], self.ph_L[ponteiro])
            fil = np.exp(lnphil) * (self.x[:,ponteiro] * Pb[ponteiro][np.newaxis,:])
            phiv = np.exp(lnphiv)
            phil = np.exp(lnphil)

            dlnphildP = self.get_dlnphidP(PR, self.x[:,ponteiro], Pb[ponteiro], self.ph_L[ponteiro])
            dlnfildP = PR.dlnfij_dP(Pb[ponteiro], dlnphildP)
            dfildP = fil * dlnfildP
            dlnphivdP = self.get_dlnphidP(PR, y[:,ponteiro], Pb[ponteiro], self.ph_V[ponteiro])
            dphivdP = phiv * dlnphivdP

            f = np.sum(fil/phiv, axis = 0) - Pb[ponteiro]
            df = np.sum((phiv * dfildP - fil * dphivdP) / phiv**2, axis=0) - 1.

            i += 1
            if i > 100 or any(Pb < 0):
                Pb[ponteiro] = 2 * self.P[ponteiro]
                ponteiro[ponteiro] = False
                break
                print("Not questionconverged - assuming its gas")

            if any(df == 0):
                import pdb; pdb.set_trace()
                raise ValueError('Change Pguess - not converging')
            #import pdb; pdb.set_trace()
            Pb[ponteiro] = Pb[ponteiro] - f / df
            K[:,ponteiro] = phil / phiv
            stop_criteria = abs(Pb[ponteiro] - Pb_old)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria <= .01*6894.757] = False
            ponteiro[ponteiro] = ponteiro_aux

        import pdb; pdb.set_trace()
        '''Y = np.copy(y)
        lnphil = np.copy(y)
        Pb = 6.8e6*np.ones(len(self.P))
        #Pb = np.sum(z* self.Pv[:,np.newaxis], axis=0)
        K = np.exp(5.37 * (1 + self.w) * (1 - 1 / (self.T / self.Tc)), dtype=np.double)[:,np.newaxis] / \
                (Pb / self.Pc[:,np.newaxis])

        while any(ponteiro):
            Y[:,ponteiro] = K[:,ponteiro] * z[:,ponteiro]
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)
            Y_old = np.copy(Y[:,ponteiro])

            Pb_old = np.copy(Pb[ponteiro])
            lnphiv = self.lnphi_based_on_deltaG(PR, y[:,ponteiro], Pb[ponteiro], self.ph_V[ponteiro])
            lnphil = self.lnphi_based_on_deltaG(PR, self.x[:,ponteiro], Pb[ponteiro], self.ph_L[ponteiro])
            fil = np.exp(lnphil) * (self.x[:,ponteiro] * Pb[ponteiro][np.newaxis,:])
            phiv = np.exp(lnphiv)
            phil = np.exp(lnphil)

            dlnphildP = self.get_dlnphidP(PR, self.x[:,ponteiro], Pb[ponteiro], self.ph_L[ponteiro])
            dlnfildP = PR.dlnfij_dP(Pb[ponteiro], dlnphildP)
            dfildP = fil * dlnfildP
            dlnphivdP = self.get_dlnphidP(PR, y[:,ponteiro], Pb[ponteiro], self.ph_V[ponteiro])
            dphivdP = phiv * dlnphivdP
            dphildP = phil * dlnphildP
            import pdb; pdb.set_trace()
            Q2 = 1 - np.sum(Y[:,ponteiro], axis=0)
            dQ2dP = - np.sum(z[:,ponteiro] * (phiv * dphildP - phil * dphivdP) / phiv ** 2, axis = 0)
            K[:,ponteiro] = phil / phiv
            Pb[ponteiro] = Pb[ponteiro] - Q2/dQ2dP
            stop_criteria = abs(Pb[ponteiro] - Pb_old)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria <= .5*6894.757] = False
            ponteiro[ponteiro] = ponteiro_aux'''

        L = self.L[ponteiro_save]
        V = self.V[ponteiro_save]
        L[self.P[ponteiro_save] > Pb[ponteiro_save]] = 1
        V[self.P[ponteiro_save] > Pb[ponteiro_save]] = 0.
        L[self.P[ponteiro_save] < Pb[ponteiro_save]] = 0.
        V[self.P[ponteiro_save] < Pb[ponteiro_save]] = 1
        self.L[ponteiro_save] = L
        self.V[ponteiro_save] = V

    def TPD(self, z): #ainda não sei onde usar isso
        x = np.zeros(self.Nc)

        #**********************Tangent Plane distance plot*********************#
        t = np.linspace(0.01, 0.99, 0.9 / 0.002) #vetor auxiliar
        TPD = np.zeros(len(t)) ##F

        for i in range(0, len(t)):
            aux = 0;
            lnphiz = self.lnphi(z, 1) #original phase

            #x = np.array([1-t[i],t[i]]) #new phase composition (1-t e t) - apenas válido para Nc=2 acredito eu.
            for k in range(0, self.Nc - 1):
                x[k] = (1 - t[i]) / (self.Nc - 1)
                x[self.Nc - 1] = t[i]

            '''O modo que x varia implica no formato de TPD. No presente exemplo,
            a fração molar do segundo componente de x varia direto com t, que é a
            variável de plotagem. Logo, a distancia dos planos tangentes será
            zero em z[Nc-1]. O contrário ocorreria'''
            lnphix = self.lnphi(x, 0); #new phase (vapor- ph=2)
            for j in range(0,self.Nc):
                fix = math.exp(lnphix[j]) * x[j] * self.P
                fiz = math.exp(lnphiz[j]) * z[j] * self.P
                aux = aux + x[j] * self.R * self.T * (math.log(fix / fiz))
                TPD[i] = aux

        plt.figure(0)
        plt.plot(t, TPD)
        plt.xlabel('x')
        plt.ylabel('TPD')
        plt.show()
        return TPD

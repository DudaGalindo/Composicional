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

    def __init__(self, w, Bin, R, Tc, Pc, T, P):
        self.w = w
        self.Bin = Bin
        self.R = R
        self.Tc = Tc
        self.Pc = Pc
        self.T = T
        self.P = P
        self.Nc = len(w)
        self.ph_L = np.ones(len(self.P), dtype = bool)
        self.ph_V = np.zeros(len(self.P), dtype = bool)
        #StabilityCheck.TPD(self)

    def run(self, z, Mw):
        PR = PengRobinson(self)
        self.equilibrium_ratio_Wilson()
        ponteiro = np.zeros(len(self.P), dtype = bool)
        dir_flash = np.argwhere(z <= 0)
        ponteiro[dir_flash[:,0]] = True
        self.L = np.empty(len(self.P))
        self.V = np.empty(len(self.P))
        self.x = np.empty(z.shape)
        self.y = np.empty(z.shape)

        if any(ponteiro):
            self.molar_properties(PR, z, ponteiro) #having a problem here
        else:
            sp1,sp2 = self.Stability(PR, z)
            ponteiro[(np.round(sp1,8) > 1) + (np.round(sp2,8) > 1)] = True
            import pdb; pdb.set_trace()
            ponteiro_save = np.copy(ponteiro)
            if any(ponteiro_save):
                self.molar_properties(PR, z, ponteiro) #having a problem here
            if any(~ponteiro_save): #tiny manipulation
                ponteiro = ~ponteiro
                ponteiro[dir_flash[:,0]] = False
                self.x[ponteiro] = 1; self.y[ponteiro] = np.copy(self.x)
                self.bubble_point_pressure()
                self.L[self.P[ponteiro] > self.Pbubble] = 1
                self.V[self.P[ponteiro] > self.Pbubble] = 0
                self.L[self.P[ponteiro] < self.Pbubble] = 0
                self.V[self.P[ponteiro] < self.Pbubble] = 1

                lnphil = self.lnphi_based_on_deltaG(PR, self.x, self.P, self.ph_L)
                lnphiv = self.lnphi_based_on_deltaG(PR, self.y, self.P, self.ph_V)
                import pdb; pdb.set_trace()
        '''self.x = z; self.y = z
        self.Mw_L, self.ksi_L, self.rho_L = self.other_properties(PR, self.x,Mw)
        self.Mw_V, self.ksi_V, self.rho_V = self.other_properties(PR, self.y,Mw)'''
        '''if sp1<1 and sp2<1:
                TPD = obj.TPD(z)
                if TPD.any()<0: #checar se isso iria funcionar
                    obj.molar_properties(z,Mw)'''

    def get_other_properties():
        PR = PengRobinson(self)
        self.Mw_L, self.ksi_L, self.rho_L = self.other_properties(PR, self.x,Mw)
        self.Mw_V, self.ksi_V, self.rho_V = self.other_properties(PR, self.y,Mw)

    def equilibrium_ratio_Wilson(self):
        self.K = np.exp(5.37 * (1 + self.w) * (1 - self.Tc / self.T)) * \
                (self.Pc  / self.P[:,np.newaxis])


    """-------------Below starts stability test calculation -----------------"""

    def Stability(self, PR, z):
        ''' In the lnphi function: 0 stands for vapor phase and 1 for liquid '''
    #****************************INITIAL GUESS******************************#
    ## Both approaches bellow should be used in case the phase is in the critical region

    #*****************************Test one**********************************#
        #Used alone when the phase investigated (z) is clearly vapor like (ph = 0)

        ponteiro = np.ones(len(self.P), dtype = bool)
        Y = z / self.K
        y = Y / np.sum(Y, axis = 1)[:, np.newaxis]
        lnphiz = PR.lnphi(self, z, self.P, self.ph_V)

        while any(ponteiro == True):
        #while max(abs(Y / Yold - 1)) > 1e-9: #convergência
            Y_old = np.copy(Y[ponteiro])
            lnphiy = PR.lnphi(self, y[ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            Y[ponteiro] = np.exp(np.log(z[ponteiro]) + lnphiz[ponteiro] - lnphiy)
            y[ponteiro] = Y[ponteiro] / np.sum(Y[ponteiro], axis = 1)[:, np.newaxis]
            stop_criteria = np.max(abs(Y[ponteiro] / Y_old - 1), axis = 1)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        stationary_point1 = np.sum(Y, axis = 1)


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
        ponteiro = np.ones(len(self.P), dtype = bool)

        Y = self.K * z
        y = Y / np.sum(Y, axis = 1)[:, np.newaxis]
        lnphiz = PR.lnphi(self, z, self.P, self.ph_L)
        while any(ponteiro):
            Y_old = np.copy(Y[ponteiro])
            lnphiy = PR.lnphi(self, y[ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            Y[ponteiro] = np.exp(np.log(z[ponteiro]) + lnphiz[ponteiro] - lnphiy)
            y[ponteiro] = Y[ponteiro] / np.sum(Y[ponteiro], axis = 1)[:, np.newaxis]
            stop_criteria = np.max(abs(Y[ponteiro] / Y_old - 1), axis = 1)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        stationary_point2 = np.sum(Y, axis = 1)

        # print(sum(Y))
        # if sum(Y) <= 1: print('estavel')#estavel2 = 1
        # else: print('instavel') #estavel2 = 0

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

        return stationary_point1,stationary_point2

    """-------------Below starts biphasic flash calculations-----------------"""
    def molar_properties(self, PR, z, ponteiro):
        self.fv = 2 * np.ones((len(self.P),self.Nc)); self.fl = np.ones((len(self.P),self.Nc)) #entrar no primeiro loop
        if self.Nc <= 2: self.molar_properties_Whitson(PR, z, ponteiro)
        else: self.molar_properties_Yinghui(PR, z, ponteiro)
        self.L = 1 - self.V

    def deltaG_molar(self, PR, l, ph):
        # "ph" is the phase originally intended for each phase.
        # l - any phase molar composition
        # Example: x - liquid(ph=1); y - vapor(ph=0)
        #   In order to verify the state of each new phase, in consequence,
        # the new phases stabilities, we have to identify the phase state
        # that makes the Gibbs free energy smaller
        lnphi = [PR.lnphi(self, l, 1 - ph), PR.lnphi(self, l, ph)]

        deltaG_molar = sum(l * (lnphi[1 - ph] - lnphi[ph]))
        if deltaG_molar<0: ph = 1 - ph
        return ph

    def deltaG_molar_vectorized(self, PR, l, P, ph):
        lnphi = np.empty([2, len(ph), len(self.w)])
        lnphi[0,:] = PR.lnphi(self, l, P, 1 - ph)
        lnphi[1,:] = PR.lnphi(self, l, P, ph)

        deltaG_molar = np.sum(l * (lnphi[1 - ph,np.arange(len(ph)),:] - lnphi[1*ph, np.arange(len(ph)),:]), axis = 1)
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

    def solve_objective_function_Yinghui(self, z1, zi, z, K1, KNc, Ki, K, x, z1_zero):
        x1_min = z1 * (1 - KNc) / (K1 - KNc)
        x1_max = (1 - KNc) / (K1 - KNc)

        vols_z_neg = np.zeros(len(K1), dtype = bool)
        vols_z_neg[np.sum(z < 0, axis = 1, dtype=bool)] = True
        theta = np.ones(z[vols_z_neg].shape)
        KNc_z_neg = KNc[vols_z_neg]
        K1_z_neg = K1[vols_z_neg]
        z1_z_neg = z1[vols_z_neg]
        K_z_neg = K[vols_z_neg]
        z_z_neg = z[vols_z_neg]

        vols_z_neg_num = np.sum(vols_z_neg*1) + 1 - np.sign(np.sum(vols_z_neg*1))
        K_z_neg_K_big1 = K_z_neg[K_z_neg > 1].reshape(vols_z_neg_num,int(len(K_z_neg[K_z_neg>1])/vols_z_neg_num))
        theta = np.ones(z[vols_z_neg].shape)
        theta[K_z_neg > 1] = ((1 - KNc_z_neg[:,np.newaxis]) / (K_z_neg_K_big1 - KNc_z_neg[:,np.newaxis])).ravel()
        aux_eq = (K_z_neg - 1) * z1_z_neg[:, np.newaxis] / (z_z_neg * (K1_z_neg[:, np.newaxis] - 1) /
                theta - (K1_z_neg[:, np.newaxis] - K_z_neg))

        #aux_eq = (K - 1) * z1 / (z * (K1 - 1) / theta - (K1 - K))
        cond = (K_z_neg[z_z_neg != 0] - 1) * z1_z_neg[:,np.newaxis] / z_z_neg[z_z_neg != 0]
        cond_aux = np.ones(cond.shape[0], dtype = bool)
        cond_aux[np.sum(cond <= 0, axis = 1, dtype=bool)] = False
        aux_eq_cond = aux_eq[cond_aux]
        vols_aux = len(cond_aux==True)

        vols_aux = np.sum(cond_aux*1) + 1 - np.sign(np.sum(cond_aux*1))
        aux_eq_cond = aux_eq_cond[aux_eq_cond >= 0].reshape(vols_aux,int(len(aux_eq_cond[aux_eq_cond >= 0])/vols_aux))
        x1_max_aux = np.copy(x1_max[vols_z_neg])
        if len(cond_aux) > 0: x1_max_aux[cond_aux] = np.min(aux_eq_cond, axis = 1)
        x1_max[vols_z_neg] = x1_max_aux

        x1_min_aux = np.copy(x1_min[vols_z_neg])
        x1_min_aux[~cond_aux] = np.max(aux_eq[~cond_aux], axis = 1)
        x1_min_aux[x1_min_aux < 0] = 0
        x1_min[vols_z_neg] = x1_min_aux

        if any(x1_min > x1_max): raise ValueError('There is no physical root')

        x1 = (x1_min + x1_max) / 2

        ponteiro = np.ones(len(x1), dtype = bool)

        while any(ponteiro):
            f = 1 + ((K1[ponteiro] - KNc[ponteiro]) / (KNc[ponteiro] - 1)) * x1[ponteiro] + np.sum(((Ki[ponteiro] - KNc[ponteiro][:, np.newaxis]) /
                (KNc[ponteiro][:, np.newaxis] - 1)) * zi[ponteiro] * (K1[ponteiro][:, np.newaxis] - 1) * x1[ponteiro][:, np.newaxis]
                / ((Ki[ponteiro] - 1) * z1[ponteiro][:, np.newaxis] + (K1[ponteiro][:,np.newaxis] - Ki[ponteiro]) *
                x1[ponteiro][:,np.newaxis]), axis = 1)
            df = ((K1[ponteiro] - KNc[ponteiro]) / (KNc[ponteiro] - 1)) + np.sum(((Ki[ponteiro] - KNc[ponteiro][:, np.newaxis]) /
                (KNc[ponteiro][:, np.newaxis] - 1)) * zi[ponteiro] * z1[ponteiro][:, np.newaxis] * (K1[ponteiro][:, np.newaxis] - 1) *
                (Ki[ponteiro] - 1) / ((Ki[ponteiro] - 1) * z1[ponteiro][:, np.newaxis] + (K1[ponteiro][:, np.newaxis] - Ki[ponteiro]) *
                x1[ponteiro][:, np.newaxis]) ** 2, axis = 1)
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


        xi = (K1[:, np.newaxis] - 1) * zi * x1[:, np.newaxis] / ((Ki - 1) * z1[:, np.newaxis] +
            (K1[:, np.newaxis] - Ki) * x1[:, np.newaxis])

        x_not_z1_zero = np.zeros(x[~z1_zero].shape)
        x_not_z1_zero[K == K1[:,np.newaxis]] = x1
        x_not_z1_zero[K == KNc[:,np.newaxis]] = 1 - np.sum(xi, axis = 1) - x1
        aux_xi = np.ones(x_not_z1_zero.shape,dtype=bool)
        aux_xi[K == K1[:,np.newaxis]] = False
        aux_xi[K == KNc[:,np.newaxis]] = False
        x_not_z1_zero[aux_xi] = xi.ravel()
        x[~z1_zero] = x_not_z1_zero
        return x

    def Yinghui_method(self, z, ponteiro):

        """ Shaping K to Nc-2 components by removing K1 and KNc and z to Nc-2
        components by removing z1 and zNc """
        K = self.K[ponteiro]
        x = self.x[ponteiro]
        K1 = np.max(K, axis=1); KNc = np.min(K, axis=1)
        z1 = z[K == K1[:,np.newaxis]]

        aux = np.ones(K.shape, dtype = bool)
        aux[K == K1[:,np.newaxis]] = False
        aux[K == KNc[:,np.newaxis]] = False
        Ki = K[aux]
        zi = z[aux]

        ''' Reshaping them into the original matricial form '''
        vols_ponteiro = np.sum(ponteiro*1) + 1 - np.sign(np.sum(ponteiro*1))
        Ki = Ki.reshape(vols_ponteiro, int(len(Ki)/vols_ponteiro))
        zi = zi.reshape(vols_ponteiro, int(len(zi)/vols_ponteiro))

        #starting x

        """ Solution """
        z1_zero = np.zeros(vols_ponteiro,dtype=bool)
        z1_zero[z1 == 0] = True

        x = self.solve_objective_function_Yinghui(z1[~z1_zero], zi[~z1_zero], z[~z1_zero],
                                K1[~z1_zero], KNc[~z1_zero], Ki[~z1_zero], K[~z1_zero], x, z1_zero)

        '''Explicit Calculation of xi'''
        #self.solve_objective_function_Yinghui_explicitly()
        z_z1_zero = z[z1_zero]
        K_z1_zero = K[z1_zero]
        K_KNc_z1_zero = K_z1_zero[K_z1_zero == KNc[z1_zero][:,np.newaxis]]

        aux_xNc = np.zeros(K_z1_zero.shape, dtype = bool); aux_x1 = np.copy(aux_xNc)
        aux_xNc[K_z1_zero == KNc[z1_zero][:,np.newaxis]] = True
        aux_x1[K_z1_zero == K1[z1_zero][:,np.newaxis]] = True
        aux_xi = aux_xNc + aux_x1
        aux_xi = ~aux_xi
        xi_z1_zero = (K1[z1_zero][:, np.newaxis] - 1) * zi[z1_zero] / (K1[z1_zero][:, np.newaxis] - Ki[z1_zero])
        x_z1_zero = np.zeros(x[z1_zero].shape)
        x_z1_zero[aux_xNc] = (K1[z1_zero] - 1) * z_z1_zero[aux_xNc] / (K1[z1_zero] - K_z1_zero[aux_xNc])
        x_z1_zero[aux_x1] = 1 - np.sum(xi_z1_zero, axis = 1) - np.sum(x_z1_zero, axis = 1)
        x_z1_zero[aux_xi] = xi_z1_zero.ravel()
        x[z1_zero] = x_z1_zero
        self.x[ponteiro] = x
        self.y = self.K * self.x


    def molar_properties_Yinghui(self, PR, z, ponteiro):
        #razao = fl/fv -> an arbitrary vector to enter in the iterative mode

        razao = np.ones(z[ponteiro].shape)/2
        ponteiro_save = np.copy(ponteiro)
        while any(ponteiro):
            self.Yinghui_method(z[ponteiro], ponteiro)
            lnphil = self.lnphi_based_on_deltaG(PR, self.x[ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            lnphiv = self.lnphi_based_on_deltaG(PR, self.y[ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            self.fl[ponteiro] = np.exp(lnphil) * (self.x[ponteiro] * self.P[ponteiro][:, np.newaxis])
            self.fv[ponteiro] = np.exp(lnphiv) * (self.y[ponteiro] * self.P[ponteiro][:, np.newaxis])
            razao[ponteiro] = np.divide(self.fl[ponteiro], self.fv[ponteiro], out = razao[ponteiro] / razao[ponteiro] * (1 + 1e-10),
                              where = self.fv[ponteiro] != 0)
            self.K[ponteiro] = razao[ponteiro] * self.K[ponteiro]
            stop_criteria = np.max(abs(razao[ponteiro] - 1), axis = 1)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        V = (z[ponteiro_save][self.x[ponteiro_save] != 0] - self.x[ponteiro_save][self.x[ponteiro_save] != 0]) / \
                          (self.y[ponteiro_save][self.x[ponteiro_save] != 0] - self.x[ponteiro_save][self.x[ponteiro_save] != 0])
        
        self.V[ponteiro_save] = V.reshape(len(ponteiro_save),int(len(V)/len(ponteiro_save)))[:,0]

    def solve_objective_function_Whitson(self, z, ponteiro):
        """ Solving for V """
        K = self.K[ponteiro]
        Vmax = 1 / (1 - np.min(self.K[ponteiro], axis = 1))
        Vmin = 1 / (1 - np.max(self.K[ponteiro], axis = 1))
        #Vmin = ((K1-KNc)*z[self.K==K1]-(1-KNc))/((1-KNc)*(K1-1))
        #proposed by Li et al for Whitson method
        V = (Vmin + Vmax) / 2
        ponteiro_save = np.copy(ponteiro)
        while any(ponteiro):
            Vold = V[ponteiro]
            f = np.sum((K[ponteiro] - 1) * z[ponteiro] / (1 + V[ponteiro][:,np.newaxis] * (K[ponteiro] - 1)), axis = 1)
            df = - np.sum((K[ponteiro] - 1) ** 2 * z[ponteiro] / (1 + V[ponteiro][:,np.newaxis] * (K[ponteiro] - 1)) ** 2, axis = 1)
            V[ponteiro] = V[ponteiro] - f / df #Newton-Raphson iterative method
            V_aux = V[ponteiro]
            V_aux[V_aux > Vmax[ponteiro]] = Vmax[ponteiro][V_aux > Vmax[ponteiro]] #(Vmax + Vold)/2
            V_aux[V_aux < Vmin[ponteiro]] = Vmin[ponteiro][V_aux < Vmin[ponteiro]] #(Vmax + Vold)/2
            V[ponteiro] = V_aux
            stop_criteria = abs(V[ponteiro] / Vold - 1)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux
            #import pdb; pdb.set_trace()
        self.V[ponteiro_save] = V[ponteiro_save]
        self.x[ponteiro_save] = z[ponteiro_save] / (1 + self.V[ponteiro_save][:,np.newaxis] * (self.K[ponteiro_save] - 1))
        self.y[ponteiro_save] = self.K[ponteiro_save] * self.x[ponteiro_save]
        # Lmax = max(self.K)/(max(self.K)-1)
        # Lmin = min(self.K)/(min(self.K)-1)
        # L = (Lmin + Lmax)/2
        # Lold = L/2
        # import pdb; pdb.set_trace()
        # while abs(L / Lold - 1) > 1e-9:
        #     Lold = L
        #     f = sum((1 - self.K) * z / (Lold + (1 - Lold) * self.K ))
        #     df = -sum((1 - self.K) ** 2 * z / (Lold + (1 - Lold) * self.K) ** 2)
        #     L = Lold - f / df #Newton-Raphson iterative method
        #     if L > Lmax: L = (Lmax + Lold)/2
        #     elif L < Lmin: L = (Lmin + Lold)/2
        # # self.L = L
        # # self.V = (1 - self.L)
        # self.x = z / (L + (1 - L) * self.K)
        # self.y = self.K * self.x

    def molar_properties_Whitson(self, PR, z, ponteiro):
        razao = np.ones(z[ponteiro].shape)/2

        while any(ponteiro):
            ponteiro_ = np.copy(ponteiro)
            self.solve_objective_function_Whitson(z, ponteiro_)
            lnphil = self.lnphi_based_on_deltaG(PR, self.x[ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            lnphiv = self.lnphi_based_on_deltaG(PR, self.y[ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            self.fv[ponteiro] = np.exp(lnphiv[ponteiro]) * (self.y[ponteiro] * self.P[ponteiro][:,np.newaxis])
            self.fl[ponteiro] = np.exp(lnphil[ponteiro]) * (self.x[ponteiro] * self.P[ponteiro][:,np.newaxis])
            razao[ponteiro] = np.divide(self.fl[ponteiro], self.fv[ponteiro], out = razao[ponteiro] / razao[ponteiro] * (1 + 1e-10),
                              where = self.fv[ponteiro] != 0)
            self.K[ponteiro] = razao[ponteiro] * self.K[ponteiro]
            stop_criteria = np.max(abs(self.fv[ponteiro] / self.fl[ponteiro] - 1), axis = 1)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

    def other_properties(self, PR, l, Mw):
        #l - any phase molar composition

        A, B = PR.coefficients_cubic_EOS_all(self, l)
        #ph = self.deltaG_molar(PR, l, 1)
        Z = PengRobinson.Z_all(B, A, 1)
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

    def bubble_point_pressure(self):
        #Isso vem de uma junção da Lei de Dalton com a Lei de Raoult
        Pv = self.K * self.P[:,np.newaxis]
        self.Pbubble = np.sum(self.x * Pv, axis = 1)

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

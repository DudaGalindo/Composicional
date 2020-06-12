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
        self.ph_L = np.ones(len(self.P), dtype = bool)
        self.ph_V = np.zeros(len(self.P), dtype = bool)

        if any(ponteiro):
            self.molar_properties(PR, z, ponteiro) #having a problem here
        else:
            sp1,sp2 = self.Stability(PR, z)
            ponteiro = ~ponteiro
            flash = np.argwhere(z <= 0)
            ponteiro[np.round(sp1,8) <= 1] = False
            ponteiro[np.round(sp2,8) <= 1] = False
            if any(ponteiro):
                self.molar_properties(PR, z, ponteiro) #having a problem here
            if any(~ponteiro): #tiny manipulation
                ponteiro = ~ponteiro
                ponteiro[dir_flash[:,0]] = False
                self.x[ponteiro] = 1; self.y[ponteiro] = np.copy(self.x)
                self.bubble_point_pressure()
                self.L[self.P[ponteiro] > self.Pbubble] = 1
                self.V[self.P[ponteiro] > self.Pbubble] = 0
                self.L[self.P[ponteiro] < self.Pbubble] = 0
                self.V[self.P[ponteiro] < self.Pbubble] = 1

                lnphil = self.lnphi_based_on_deltaG(PR, self.x, self.ph_L)
                lnphiv = self.lnphi_based_on_deltaG(PR, self.y, self.ph_V)
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
        Yold = 0.9 * Y
        y = Y / np.sum(Y, axis = 1)[:, np.newaxis]
        lnphiz = PR.lnphi(self, z, self.ph_V)

        while any(ponteiro == True):
        #while max(abs(Y / Yold - 1)) > 1e-9: #convergência
            Yold = np.copy(Y[ponteiro])
            lnphiy = PR.lnphi(self, y, self.ph_L[ponteiro])
            Y = np.exp(np.log(z[ponteiro]) + lnphiz[ponteiro] - lnphiy)
            y = Y / np.sum(Y, axis = 1)[:, np.newaxis]
            stop_criteria = np.max(abs(Y / Yold - 1), axis = 1)
            ponteiro = ponteiro[ponteiro]
            ponteiro[stop_criteria < 1e-9] = False
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
        Y_old = 0.9 * Y
        y = Y / np.sum(Y, axis = 1)[:, np.newaxis]
        lnphiz = PR.lnphi(self, z, np.ones(len(self.P)))
        while any(ponteiro):
            Y_old = np.copy(Y)
            lnphiy = PR.lnphi(self, y.T, np.zeros(len(self.P))[ponteiro])
            Y = np.exp(np.log(z[ponteiro]) + lnphiz[ponteiro] - lnphiy)
            y = Y / np.sum(Y, axis = 1)[:, np.newaxis]
            stop_criteria = np.max(abs(Y / Yold - 1), axis = 1)
            ponteiro = ponteiro[ponteiro]
            ponteiro[stop_criteria < 1e-9] = False
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

    def deltaG_molar_vectorized(self, PR, l, ph):
        lnphi = np.empty([2, len(ph), len(self.w)])
        lnphi[0,:] = PR.lnphi(self, l, 1 - ph)
        lnphi[1,:] = PR.lnphi(self, l, ph)

        deltaG_molar = np.sum(l * (lnphi[1 - ph,np.arange(len(self.P)),:] - lnphi[1*ph, np.arange(len(self.P)),:]), axis = 1)
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


    def lnphi_based_on_deltaG(self, PR, l, ph):
        ph = self.deltaG_molar_vectorized(PR, l, ph)
        return PR.lnphi(self, l, ph)

    def solve_objective_function_Yinghui(self, z1, zi, z, K1, KNc, Ki, z1_zero):
        K = self.K[~z1_zero]
        x1_min = z1 * (1 - KNc) / (K1 - KNc)
        x1_max = (1 - KNc) / (K1 - KNc)

        vols_z_neg = np.zeros(len(K1), dtype = bool)
        vols_z_neg[np.sum(z < 0, axis = 1, dtype=bool)] = True
        #theta = np.ones(z[vols_z_neg].shape)
        import pdb; pdb.set_trace()

        if any(z < 0):
            theta = np.ones(len(z))
            theta[self.K > 1] = (1 - KNc) / (self.K[self.K > 1] - KNc)
            aux_eq = (self.K - 1) * z1 / (z * (K1 - 1) / theta - (K1 - self.K))
            if all((self.K[z != 0] - 1) * z1 / z[z != 0] > 0):
                aux_eq = aux_eq[aux_eq >= 0] #se arr<0 é descartado
                x1_max = min(aux_eq)
            else: x1_min = max(np.append(aux_eq, 0))


        if any(x1_min > x1_max): raise ValueError('There is no physical root')

        x1 = (x1_min + x1_max) / 2

        f_aux = np.ones(len(x1), dtype = bool)
        while any(f_aux):
            f = 1 + ((K1 - KNc) / (KNc - 1)) * x1 + sum(((Ki - KNc) / (KNc - 1))
               * zi * (K1 - 1) * x1 / ((Ki - 1) * z1 + (K1[:,np.newaxis] - Ki) * x1[:,np.newaxis]))
            df = ((K1 - KNc) / (KNc - 1)) + sum(((Ki - KNc) / (KNc - 1)) * zi *
                z1 * (K1 - 1)* (Ki - 1) / ((Ki - 1) * z1 + (K1 - Ki) * x1) ** 2)
            x1 = x1 - f/df #Newton-Raphson iterative method
            x1[x1 > x1_max] = (x1_min + x1_max)/2
            x1[x1 < x1_min] = (x1_min + x1_max)/2
            x1_max[f * df > 0] = x1
            x1_min[f * df < 0] = x1
            f_aux[f < 1e-10] = False


        xi = (K1 - 1) * zi * x1 / ((Ki - 1) * z1 + (K1 - Ki) * x1)
        self.x[not z1_zero][K == K1] = x1
        self.x[not z1_zero][K == KNc] = 1 - np.sum(xi, axis = 1) - x1
        return xi

    def Yinghui_method(self, z):

        """ Shaping K to Nc-2 components by removing K1 and KNc and z to Nc-2
        components by removing z1 and zNc """
        K1 = np.max(self.K, axis=1); KNc = np.min(self.K, axis=1)
        z1 = z[self.K == K1[:,np.newaxis]]

        aux = np.ones(self.K.shape, dtype = bool)
        aux[self.K == K1[:,np.newaxis]] = False
        aux[self.K == KNc[:,np.newaxis]] = False
        Ki = self.K[aux]
        zi = z[aux]

        ''' Reshaping them into the original matricial form '''
        Ki = Ki.reshape(len(self.P), int(len(Ki)/len(self.P)))
        zi = zi.reshape(len(self.P), int(len(zi)/len(self.P)))

        #starting x
        self.x = np.zeros(z.shape)
        xi = np.empty(zi.shape)
        """ Solution """
        z1_zero = np.zeros(len(self.P), dtype=bool)
        z1_zero[z1 == 0] = True

        xi[~z1_zero] = self.solve_objective_function_Yinghui(z1[~z1_zero], zi[~z1_zero],
                        z[~z1_zero], K1[~z1_zero], KNc[~z1_zero], Ki[~z1_zero], z1_zero)
        import pdb; pdb.set_trace()
        '''Explicit Calculation of xi'''
        xi[z1_zero] = (K1[z1_zero] - 1) * zi[z1_zero] / (K1[z1_zero] - Ki[z1_zero])
        self.x[self.K == KNc[:,np.newaxis]][z1_zero] = (K1[z1_zero] - 1) * z[self.K == KNc[:,np.newaxis]][z1_zero] / (K1[z1_zero] -
                                self.K[self.K == KNc[:,np.newaxis]][z1_zero])
        self.x[self.K == K1[:,np.newaxis]][z1_zero] = 1 - np.sum(xi[z1_zero], axis = 1) - np.sum(self.x[z1_zero], axis = 1)

        #ainda não sei como tirar esse for
        for i in range(0, len(xi[0,:])):
            self.x[self.K == Ki[:,i]] = xi[:,i]

        self.y = self.K * self.x

    def molar_properties_Yinghui(self, PR, z, ponteiro):
        #razao = fl/fv -> an arbitrary vector to enter in the iterative mode
        razao = np.ones(self.Nc)/2
        while any(ponteiro):
            self.Yinghui_method(z[ponteiro])
            lnphil = self.lnphi_based_on_deltaG(PR, self.x[ponteiro], self.ph_L[ponteiro])
            lnphiv = self.lnphi_based_on_deltaG(PR, self.y[ponteiro], self.ph_V[ponteiro])
            self.fl = np.exp(lnphil) * (self.x[ponteiro] * self.P[ponteiro])
            self.fv = np.exp(lnphiv) * (self.y[ponteiro] * self.P[ponteiro])
            razao = np.divide(self.fl, self.fv, out = razao / razao * (1 + 1e-10),
                              where = self.fv != 0)
            import pdb; pdb.set_trace()
            self.K = razao * self.K[ponteiro]
            stop_criteria = np.max(abs(self.fv / self.fl - 1), axis=0)
            stop = np.argwhere(stop_criteria > 1e-9)
            ponteiro[ponteiro][stop] = False

        self.V = (z[self.x != 0] - self.x[self.x != 0]) / (self.y[self.x != 0] -
                self.x[self.x != 0])
        self.V = self.V[0]

    def solve_objective_function_Whitson(self, z, ponteiro):
        """ Solving for V """
        Vmax = 1 / (1 - np.min(self.K, axis = 1))
        Vmin = 1 / (1 - np.max(self.K, axis = 1))
        #Vmin = ((K1-KNc)*z[self.K==K1]-(1-KNc))/((1-KNc)*(K1-1))
        #proposed by Li et al for Whitson method
        self.V = (Vmin + Vmax) / 2
        Vold = self.V / 2 #just to get into the loop

        while any(ponteiro):
            Vold = self.V
            f = np.sum((self.K - 1) * z / (1 + self.V * (self.K - 1)), axis = 1)
            df = -np.sum((self.K - 1) ** 2 * z / (1 + self.V * (self.K - 1)) ** 2, axis = 1)
            self.V = self.V - f / df #Newton-Raphson iterative method
            self.V[self.V > Vmax] = Vmax #(Vmax + Vold)/2
            self.V[self.V < Vmin] = Vmin #(Vmax + Vold)/2
            stop_criteria = abs(self.V / Vold - 1)
            ponteiro[stop_criteria < 1e-8] = False

        self.x = z / (1 + self.V * (self.K - 1))
        self.y = self.K * self.x
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
        razao = np.ones(self.Nc)/2
        self.ph_L = np.ones(len(self.P))
        self.ph_V = np.zeros(len(self.P))
        import pdb; pdb.set_trace()
        while any(ponteiro):
            self.solve_objective_function_Whitson(z[ponteiro])
            lnphil = self.lnphi_based_on_deltaG(PR, self.x[ponteiro], self.ph_L[ponteiro])
            lnphiv = self.lnphi_based_on_deltaG(PR, self.y[ponteiro], self.ph_V[ponteiro])
            self.fv = np.exp(lnphiv) * (self.y[ponteiro] * self.P[ponteiro])
            self.fl = np.exp(lnphil) * (self.x[ponteiro] * self.P[ponteiro])
            razao = np.divide(self.fl, self.fv, out = razao / razao * (1 + 1e-10),
                              where = self.fv != 0)
            self.K = razao * self.K[ponteiro]
            stop_criteria = np.max(abs(self.fv / self.fl - 1), axis=0)
            stop = np.argwhere(stop_criteria > 1e-9)
            ponteiro[ponteiro][stop] = False

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

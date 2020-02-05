"""Check stability of a thermodynamic equilibrium."""
import numpy as np
import math
from scipy.misc import derivative
from sympy import *
import sympy
from mpmath import *
import matplotlib.pyplot as plt

## Encontrar os pontos estacionarios. Estes correspondem aos pontos nos quais a derivada de g com respeito a Y é 0
## Todas as equações foram confirmadas pelo livro do Dandekar e a biblioteca do thermo

class StabilityCheck:
    """ Check for stability of a thermodynamic equilibrium and returns the
    equilibrium phase compositions (perform the flash calculation). """

    def __init__(self, w, Bin, R, Tc, Pc, T, P, C7):
        self.w = w
        self.Bin = Bin
        self.R = R
        self.Tc = Tc
        self.Pc = Pc
        self.T = T
        self.P = P
        self.Nc = len(w)
        self.C7 = C7
        self.n_phases = 2 #number of phases
        #StabilityCheck.TPD(self)

    def d_dn_all_derivatives(self, Nphase):
        l = np.zeros([self.Nc, self.n_phases])
        l[:,0] = self.y; l[:,1] = self.x # if the number of phase changes, this would change as well
        nkphase = l * Nphase
        delta = 0.0001
        dlnphi_dn = np.zeros([self.Nc, self.Nc, self.n_phases])
        # matrix = np.copy(dlnphi_dn)
        dZ_dn = np.zeros([self.Nc, self.n_phases])
        for ph in range(self.n_phases):
            if Nphase[ph] != 0:
                for i in range(0,self.Nc):
                    l_plus = np.copy(l[:,ph]); l_minus = np.copy(l[:,ph])
                    l_plus[i] = (nkphase[i,ph] + delta/2)/Nphase[ph]
                    l_minus[i] = (nkphase[i,ph] - delta/2)/Nphase[ph]
                    dlnphi_dn[:,i,ph] = (self.lnphi(l_plus, ph)
                     - self.lnphi(l_minus, ph))/delta
                    dZ_dn[i,ph] = (self.Z(l_plus, ph)
                     - self.Z(l_minus, ph))/delta
        return dlnphi_dn, dZ_dn

    def d_dP_all_derivatives(self, Nphase):
        l = np.zeros([self.Nc, self.n_phases])
        l[:,0] = self.y; l[:,1] = self.x
        delta = 0.0001

        P = np.copy(self.P); P_plus = np.copy(P); P_minus = np.copy(P)
        dZ_dP = np.zeros(self.n_phases)
        dlnphi_dP = np.zeros([self.Nc, self.n_phases])
        for ph in range(0, self.n_phases):
            if Nphase[ph] != 0:
                self.P = P + delta/2
                Z_plus = self.Z(l[:,ph], ph)
                lnphi_plus = self.lnphi(l[:,ph], ph)
                self.P = P - delta/2
                dZ_dP[ph] = (Z_plus - self.Z(l[:,ph], ph))/delta
                dlnphi_dP[:,ph] = (lnphi_plus - self.lnphi(l[:,ph], ph)) / delta
        self.P = P #voltando ao original
        return dlnphi_dP, dZ_dP

    def Z(self, l, ph):
        A, B = self.coefficientsPR(l)
        return StabilityCheck.Z_PR(B, A, ph)

    def dVt_derivatives(self, Nphase_allvolumes, eta_j):
        n_vols = len(Nphase_allvolumes[0,0,:])
        """ Initializing dN derivatives """
        dnij_dNk = np.zeros([self.Nc, self.Nc, self.n_phases, n_vols])
        dZj_dnij = np.zeros([self.Nc, self.n_phases, n_vols])

        """ Initializing dP derivatives """
        dnij_dP = np.zeros([self.Nc, self.n_phases, n_vols])
        dZj_dP = np.zeros([1, self.n_phases, n_vols])
        Zj = np.zeros([1,self.n_phases, n_vols])

        for b in range(n_vols):
            Nphase = Nphase_allvolumes[0,:,b]
            dlnphiij_dnij, dZj_dnij[:,:,b] = self.d_dn_all_derivatives(Nphase)
            dlnphij_dP, dZj_dP[0,:,b] = self.d_dP_all_derivatives(Nphase)
            Zj[0,:,b] = np.array([self.Z(self.y,0), self.Z(self.x,1)])
            matrix = np.sum(dlnphiij_dnij, axis = 2)

            dlnphi_dP_vector = dlnphij_dP[:,0] - dlnphij_dP[:,1]
            dnij_dP[:,1,b] =  np.linalg.inv(matrix)@dlnphi_dP_vector
            dnij_dP[:,0,b] = - dnij_dP[:,1,b]
            for k in range(self.Nc):
                dlnphil_dnk = dlnphiij_dnij[:,k,1]
                dnij_dNk[:,k,1,b] = np.linalg.inv(matrix)@dlnphil_dnk
                dnij_dNk[:,k,0,b] = - dnij_dNk[:,k,1,b]
                dnij_dNk[k,k,0,b] = 1 - dnij_dNk[k,k,1,b]

        dvj_dnij = self.R * self.T / self.P * dZj_dnij
        dVj_dNk = np.sum(dnij_dNk * (1/eta_j + Nphase[:,np.newaxis] * dvj_dnij)[np.newaxis,:,:,:], axis = 0)
        dVt_dNk = np.sum(dVj_dNk, axis = 1)

        dvj_dP = self.R * self.T / self.P * (dZj_dP - Zj / self.P)
        dVj_dP = np.sum(Nphase_allvolumes * dvj_dP + np.sum(dnij_dP * (1/eta_j + Nphase[:,np.newaxis] * dvj_dnij), axis = 0),axis = 0)
        dVt_dP = np.sum(dVj_dP,axis = 0)
        return dVt_dNk, dVt_dP

    def run(self, z, Mw):

        self.equilibrium_ratio_Wilson()
        if any(z <= 0):
            self.molar_properties(z)
        else:
            sp1,sp2 = self.Stability(z)
            if sp1 > 1 or sp2 > 1:
                self.molar_properties(z)
            else: #tiny manipulation
                self.x = np.ones(self.Nc); self.y = np.copy(self.x)
                self.bubble_point_pressure()
                if self.P > self.Pbubble: self.L = 1; self.V = 0
                else: self.L = 0; self.V = 1
                lnphil = self.lnphi_based_on_deltaG(self.x, 1)
                lnphiv = self.lnphi_based_on_deltaG(self.y, 0)
                self.fv = np.exp(lnphiv) * (self.y * self.P)
                self.fl = np.exp(lnphil) * (self.x * self.P)
                self.K = self.y/self.x
                self.z = self.x*self.L + self.y*self.V

        self.Mw_L, self.eta_L, self.rho_L = self.other_properties(self.x,Mw)
        self.Mw_V, self.eta_V, self.rho_V = self.other_properties(self.y,Mw)

        V = np.zeros([1,2,99])
        V[:,1,:] = np.ones(99)*11.33
        eta = np.zeros([1,2,99])
        eta[:,0,:] = self.eta_V; eta[:,1,:] = self.eta_L
        self.Nphase = eta * V
        self.dVt_derivatives(self.Nphase, eta)
        # import pdb; pdb.set_trace()
        # self.run_sym(z,1)
        '''if sp1<1 and sp2<1:
                TPD = obj.TPD(z)
                if TPD.any()<0: #checar se isso iria funcionar
                    obj.molar_properties(z,Mw)'''
    def Z_PR_sym(B, A, ph):
        Z = symbols('Z')
        a = Z**3 - (1-B)*Z**2 + (A-2*B-3*B**2)*Z-(A*B-B**2-B**3)
        coef = [1, -(1 - B), (A - 2*B - 3*B**2), -(A*B - B**2 - B**3)]
        Za = solve(a,Z)
        return Za

    def lnphi_sym(self, nkphase, Nphase, ph):
        if Nphase == 0:
            lnphi = np.zeros(self.Nc)
        else:
            l = nkphase/Nphase
            A, B = self.coefficientsPR(l)
            Z = StabilityCheck.Z_PR_sym(B, A, ph)

            lnphi1 =  self.b / self.bm * (Z[0] - 1) - sympy.log(Z[0] - B) - A / (2 * (2 ** (1/2))
                        * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * sympy.log((Z[0] + (1 +
                        2 ** (1/2)) * B) / (Z[0] + (1 - 2 ** (1/2)) * B))
            lnphi2 = self.b / self.bm * (Z[1] - 1) - sympy.log(Z[1] - B) - A / (2 * (2 ** (1/2))
                        * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * sympy.log((Z[1] + (1 +
                        2 ** (1/2)) * B) / (Z[1] + (1 - 2 ** (1/2)) * B))
            lnphi3 =  self.b / self.bm * (Z[2] - 1) - sympy.log(Z[2] - B) - A / (2 * (2 ** (1/2))
                        * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * sympy.log((Z[2] + (1 +
                        2 ** (1/2)) * B) / (Z[2] + (1 - 2 ** (1/2)) * B))
            lnphi = np.array([lnphi1,lnphi2,lnphi3])

        return lnphi

    def run_sym(self, z, ph):
        nkphase = sympy.symarray('nkphase', len(z))
        Nphaseaa = self.Nphase[0,:,0]
        Nphase = Nphaseaa[ph]
        if Nphase!=0:
            lnphi_all = self.lnphi_sym(nkphase, Nphase, ph)
            A, B = self.coefficientsPR(z)
            Z = np.array(StabilityCheck.Z_PR_sym(B, A, ph))
            reais = np.full(len(Z), False, dtype=bool)
            for r in range(0,len(Z)):
                reais[r] = sympify(Z[r]).is_real

            Z_reais = Z[reais]
            Z_ans = min(Z_reais) * ph + max(Z_reais) * (1 - ph)
            a = np.argwhere(Z == Z_ans)

            lnphi = lnphi_all[a,:]
            lnphi = sympy.simplify(lnphi)
            dlnphi_dnk = np.zeros([self.Nc,self.Nc])
            nkphase_value = z * Nphase
            for i in range(0,self.Nc):
                dlnphi_dn = sympy.diff(lnphi,nkphase[i])
                import pdb; pdb.set_trace()
                func1 = lambdify(nkphase, dlnphi_dn,'numpy')
                dlnphi_dnk[:,i] = np.array(func1(*nkphase_value)[0][0])
        else: dlnphi_dnk = np.zeros([self.Nc,self.Nc])

        return dlnphi_dnk

    def equilibrium_ratio_Wilson(self):
        self.K = np.exp(5.37 * (1 + self.w) * (1 - self.Tc / self.T)) * \
                (self.Pc / self.P)

    def coefficientsPR(self, l):
        #l - any phase molar composition
        PR_kC7 = np.array([0.379642, 1.48503, 0.1644, 0.016667])
        PR_k = np.array([0.37464, 1.5422, 0.26992])

        k = (PR_kC7[0] + PR_kC7[1] * self.w - PR_kC7[2] * self.w ** 2 + \
            PR_kC7[3] * self.w ** 3) * self.C7 + (PR_k[0] + PR_k[1] * self.w - \
            PR_k[2] * self.w ** 2) * (1 - self.C7)
        alpha = (1 + k * (1 - (self.T / self.Tc) ** (1 / 2))) ** 2
        aalpha_i = 0.45724 * (self.R * self.Tc) ** 2 / self.Pc * alpha
        self.b = 0.07780 * self.R * self.Tc / self.Pc
        aalpha_i_reshape = np.ones((self.Nc,self.Nc)) * aalpha_i[:,np.newaxis]
        aalpha_ij = np.sqrt(aalpha_i_reshape.T * aalpha_i[:,np.newaxis]) \
                        * (1 - self.Bin)
        self.bm = sum(l * self.b)
        B = self.bm * self.P / (self.R * self.T)
        l_reshape = np.ones((aalpha_ij).shape) * l[:, np.newaxis]
        self.aalpha = (l_reshape.T * l[:,np.newaxis] * aalpha_ij).sum()
        A = self.aalpha * self.P / (self.R * self.T) ** 2
        self.psi = (l_reshape * aalpha_ij).sum(axis = 0)

        return A, B

    def Z_PR(B, A, ph):
        # PR cubic EOS: Z**3 - (1-B)*Z**2 + (A-2*B-3*B**2)*Z-(A*B-B**2-B**3)
        coef = [1, -(1 - B), (A - 2*B - 3*B**2), -(A*B - B**2 - B**3)]
        Z = np.roots(coef)
        root = np.isreal(Z) # return True for real roots
        #position where the real roots are - crated for organization only
        real_roots_position = np.where(root == True)
        Z_reais = np.real(Z[real_roots_position[:]]) #Saving the real roots
        Z_ans = min(Z_reais) * ph + max(Z_reais) * (1 - ph)

        ''' This last line, considers that the phase is composed by a pure
         component, so the EOS model can return more than one real root.
            If liquid, Zl = min(Z) and gas, Zv = max(Z).
            You can notice that, if there's only one real root,
        it works as well.'''
        return Z_ans

    def lnphi(self, l, ph):
        #l - any phase molar composition
        A, B = self.coefficientsPR(l)
        Z = StabilityCheck.Z_PR(B, A, ph)
        lnphi = (self.b / self.bm )* (Z - 1) - np.log(Z - B) - A / (2 * (2 ** (1/2))
                * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * np.log((Z + (1 +
                2 ** (1/2)) * B) / (Z + (1 - 2 ** (1/2)) * B))

        return lnphi

    """-------------Below starts stability test calculation -----------------"""

    def Stability(self, z):
        ''' In the lnphi function: 0 stands for vapor phase and 1 for liquid '''
    #****************************INITIAL GUESS******************************#
    ## Both approaches bellow should be used in case the phase is in the critical region

    #*****************************Test one**********************************#
        #Used alone when the phase investigated (z) is clearly vapor like (ph = 0)

        Y = z / self.K
        Yold = 0.9 * Y
        y = Y / sum(Y)
        lnphiz = self.lnphi(z, 0)
        while max(abs(Y / Yold - 1)) > 1e-9: #convergência
            Yold = np.copy(Y)
            lnphiy = self.lnphi(y, 1)
            Y = np.exp(np.log(z) + lnphiz - lnphiy)
            y = Y / sum(Y);

        stationary_point1 = sum(Y)
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

        Y = self.K * z
        Y_old = 0.9 * Y
        y = Y / sum(Y)
        lnphiz = self.lnphi(z, 1)
        while max(abs(Y / Y_old - 1)) > 1e-9:
            Y_old = np.copy(Y)
            lnphiy = self.lnphi(y, 0)
            Y = np.exp(np.log(z) + lnphiz - lnphiy)
            y = Y / sum(Y)

        stationary_point2 = sum(Y)
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
    def molar_properties(self, z):
        self.fv = 2 * np.ones(self.Nc); self.fl = np.ones(self.Nc) #entrar no primeiro loop
        if self.Nc <= 2: self.molar_properties_Whitson(z)
        else: self.molar_properties_Yinghui(z)
        self.L = 1 - self.V

    def deltaG_molar(self, l, ph):
        # "ph" is the phase originally intended for each phase.
        # l - any phase molar composition
        # Example: x - liquid(ph=1); y - vapor(ph=0)
        #   In order to verify the state of each new phase, in consequence,
        # the new phases stabilities, we have to identify the phase state
        # that makes the Gibbs free energy smaller

        lnphi = [self.lnphi(l, 1 - ph), self.lnphi(l, ph)]
        deltaG_molar = sum(l * (lnphi[1 - ph] - lnphi[ph]))
        if deltaG_molar >= 0: ph = ph
        else: ph = 1 - ph
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

    def lnphi_based_on_deltaG(self, l, ph):
        ph = self.deltaG_molar(l, ph)
        return self.lnphi(l,ph)

    def solve_objective_function_Yinghui(self, z1, zi, z, K1, KNc, Ki):

        x1_min = z1 * (1 - KNc) / (K1 - KNc)
        x1_max = (1 - KNc) / (K1 - KNc)

        if any(z < 0):
            theta = np.ones(len(z))
            theta[self.K > 1] = (1 - KNc) / (self.K[self.K > 1] - KNc)
            aux_eq = (self.K - 1) * z1 / (z * (K1 - 1) / theta - (K1 - self.K))
            if all((self.K[z != 0] - 1) * z1 / z[z != 0] > 0):
                aux_eq = aux_eq[aux_eq >= 0] #se arr<0 é descartado
                x1_max = min(aux_eq)
            else: x1_min = max(np.append(arr, 0))

        if x1_min > x1_max: raise ValueError('There is no physical root')

        x1 = (x1_min + x1_max) / 2
        f = 1
        while abs(f) > 1e-10:
            f = 1 + ((K1 - KNc) / (KNc - 1)) * x1 + sum(((Ki - KNc) / (KNc - 1))
               * zi * (K1 - 1) * x1 / ((Ki - 1) * z1 + (K1 - Ki) * x1))
            df = ((K1 - KNc) / (KNc - 1)) + sum(((Ki - KNc) / (KNc - 1)) * zi *
                z1 * (K1 - 1)* (Ki - 1) / ((Ki - 1) * z1 + (K1 - Ki) * x1) ** 2)
            x1 = x1 - f/df #Newton-Raphson iterative method
            if (x1 > x1_max) or (x1 < x1_min): x1 = (x1_min + x1_max)/2 #recommended
            if f * df > 0: x1_max = x1
            if f * df < 0: x1_min = x1

        xi = (K1 - 1) * zi * x1 / ((Ki - 1) * z1 + (K1 - Ki) * x1)
        self.x[self.K == K1] = x1
        self.x[self.K == KNc] = 1 - sum(xi) - x1
        return xi

    def Yinghui_method(self, z):

        """ Shaping K to Nc-2 components by removing K1 and KNc """
        K1 = max(self.K); KNc = min(self.K)
        Ki = self.K[(self.K != K1) & (self.K != KNc)]

        """ Shaping z to Nc-2 components by removing z1 and zNc """
        z1 = z[self.K == K1]
        index1 = np.argwhere(self.K == K1)
        zi = np.delete(z, index1)
        indexNc = np.argwhere(self.K == KNc)
        zi = np.delete(zi, indexNc - 1)

        #starting x
        self.x = np.zeros(self.Nc)

        """ Solution """
        if z1 != 0:
            xi = self.solve_objective_function_Yinghui(z1, zi, z, K1, KNc, Ki)
        else:
            '''Explicit Calculation of xi'''
            xi = (K1 - 1) * zi / (K1 - Ki)
            self.x[self.K == KNc] = (K1 - 1) * z[self.K == KNc] / (K1 -
                                    self.K[self.K == KNc])
            self.x[self.K == K1] = 1 - sum(xi) - sum(self.x)

        #ainda não sei como tirar esse for
        for i in range(0, len(xi)):
            self.x[self.K == Ki[i]] = xi[i]

        self.y = self.K * self.x

    def molar_properties_Yinghui(self, z):
        #razao = fl/fv -> an arbitrary vector to enter in the iterative mode
        razao = np.ones(self.Nc)/2
        while max(abs(razao - 1)) > 1e-9:
            self.Yinghui_method(z)
            lnphil = self.lnphi_based_on_deltaG(self.x, 1)
            lnphiv = self.lnphi_based_on_deltaG(self.y, 0)
            self.fl = np.exp(lnphil) * (self.x * self.P)
            self.fv = np.exp(lnphiv) * (self.y * self.P)
            razao = np.divide(self.fl, self.fv, out = razao / razao * (1 + 1e-10),
                              where = self.fv != 0)
            self.K = razao * self.K

        self.V = (z[self.x != 0] - self.x[self.x != 0]) / (self.y[self.x != 0] -
                self.x[self.x != 0])
        self.V = self.V[0]

    def solve_objective_function_Whitson(self, z):
        """ Solving for V """
        Vmax = 1 / (1 - min(self.K))
        Vmin = 1 / (1 - max(self.K))
        #Vmin = ((K1-KNc)*z[self.K==K1]-(1-KNc))/((1-KNc)*(K1-1))
        #proposed by Li et al for Whitson method
        self.V = (Vmin + Vmax) / 2
        Vold = self.V / 2 #just to get into the loop

        while abs(self.V / Vold - 1) > 1e-8:
            Vold = self.V
            f = sum((self.K - 1) * z / (1 + self.V * (self.K - 1)))
            df = -sum((self.K - 1) ** 2 * z / (1 + self.V * (self.K - 1)) ** 2)
            self.V = self.V - f / df #Newton-Raphson iterative method

            if self.V > Vmax: self.V = Vmax #(Vmax + Vold)/2
            elif self.V < Vmin: self.V = Vmin #(Vmin + Vold)/2

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

    def molar_properties_Whitson(self, z):
        razao = np.ones(self.Nc)/2
        while max(abs(razao - 1)) > 1e-9:
            self.solve_objective_function_Whitson(z)
            lnphil = self.lnphi_based_on_deltaG(self.x, 1)
            lnphiv = self.lnphi_based_on_deltaG(self.y, 0)
            self.fv = np.exp(lnphiv) * (self.y * self.P)
            self.fl = np.exp(lnphil) * (self.x * self.P)
            razao = np.divide(self.fl, self.fv, out = razao / razao * (1 + 1e-10),
                              where = self.fv != 0)
            self.K = razao * self.K

    def other_properties(self, l, Mw):
        #l - any phase molar composition
        A, B = self.coefficientsPR(l)
        ph = self.deltaG_molar(l, 1)
        Z = StabilityCheck.Z_PR(B, A, ph)
        eta_phase = self.P / (Z * self.R * self.T)
        Mw_phase = sum(l * Mw)
        rho_phase = eta_phase * sum(l * Mw)
        """
        Legenda:
            rho: phase mass density
            Mw: phase molecular weight
            eta: phase molar density
        """
        # se precisar retornar mais coisa, entra aqui
        return Mw_phase, eta_phase, rho_phase

    def bubble_point_pressure(self):
        #Isso vem de uma junção da Lei de Dalton com a Lei de Raoult
        Pv = self.K*self.P
        self.Pbubble = sum(self.x * Pv)

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

"""Check stability of a thermodynamic equilibrium."""
import numpy as np
import math
#import thermo
from scipy.misc import derivative
from equation_of_state import PengRobinson
import matplotlib.pyplot as plt
import time
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
        self.pressure_water_saturation()
        self.equilibrium_ratio_aqueous(z)
        self.L = np.empty(len(self.P))
        self.V = np.empty(len(self.P))
        self.A = np.empty(len(self.P))     # 3 phases
        self.x = np.empty(z.shape)
        self.y = np.empty(z.shape)
        self.a = np.empty(z.shape)         # 3 phases

        ponteiro_flash = np.zeros(len(self.P), dtype = bool)
        dir_flash = np.argwhere(z <= 0)
        ponteiro_flash[dir_flash[:,1]] = True

        #sp = np.zeros(5) # In the case of 5 trial phases in stability analysis

        if any(~ponteiro_flash):
            sp, Kvalue = self.Stability_2phase(PR, z, np.copy(~ponteiro_flash))
            sp = np.round(sp, 14)
            #import pdb; pdb.set_trace()
            ponteiro_aux = ponteiro_flash[~ponteiro_flash]
            ponteiro_aux[(sp>1).sum(axis=0,dtype=bool)] = True #os que devem passar para o calculo de flash
            ponteiro_flash[~ponteiro_flash] = ponteiro_aux
            #index_spmax = np.argmax(np.round(sp, 10))
            #if sp[4] > 1:
                    #self.K = self.Kw
            #import pdb; pdb.set_trace()
            self.K[:,sp[4]>1] = self.Kw[:,sp[4]>1]
            #self.equilibrium_ratio_2flash(Kvalue[index_spmax])
            #print(f'sp: {sp}')
            #print(Kvalue)
            #self.K = Kvalue[1].copy()
            #print(self.K)



            #sp1,sp2 = self.Stability(PR, z, np.copy(~ponteiro_flash))
            #import pdb; pdb.set_trace()

            #ponteiro_aux = ponteiro_flash[~ponteiro_flash]
            #ponteiro_aux[(np.round(sp1,14) > 1) + (np.round(sp2,14) > 1)] = True #os que devem passar para o calculo de flash
            #ponteiro_flash[~ponteiro_flash] = ponteiro_aux
            #import pdb; pdb.set_trace()

        #import pdb; pdb.set_trace()
        self.molar_properties(PR, z, np.ones_like(ponteiro_flash, dtype=bool)) # cálculo do flash bifásico
        #import pdb; pdb.set_trace()

        self.y[:,(self.L>1) + (self.V>1)] = z[:,(self.L>1) + (self.V>1)]
        self.x[:,(self.L>1) + (self.V>1)] = z[:,(self.L>1) + (self.V>1)]
        self.L[self.L>1] = 1
        self.L[self.L<0] = 0
        self.V = 1 - self.L
        #ponteiro_flash[self.L<0] = False
        #ponteiro_flash[self.L>1] = False

        #self.x[:,~ponteiro_flash] = z[:,~ponteiro_flash]
        #self.y[:,~ponteiro_flash] = z[:,~ponteiro_flash]
        #self.bubble_point_pressure(PR, z, Mw, np.copy(~ponteiro_flash))
        #self.x = z
        #self.y = z
        #self.L = 1
        #self.V = 0
        self.get_other_properties(PR, Mw)

        '-----------------REVER----------------'
        ponteiro_flash_3phase = np.zeros(len(self.P), dtype = bool)
        ponteiro_flash_3phase[(self.L != 1) & (self.V != 1)] = True
        ponteiro_flash_3phase2 = ponteiro_flash_3phase.copy()

        sp2, Kvalue2 = self.Stability_3phase(PR, self.x, np.copy(ponteiro_flash_3phase))
        sp2 = np.round(sp2, 8)
        #import pdb; pdb.set_trace()
        ponteiro_aux = ~ponteiro_flash_3phase[ponteiro_flash_3phase]
        #ponteiro_aux = np.zeros(len(self.P), dtype = bool)
        ponteiro_aux[(sp2>1).sum(axis=0,dtype=bool)] = True
        ponteiro_flash_3phase[ponteiro_flash_3phase] = ponteiro_aux
        #ponteiro_flash_3phase = ponteiro_aux

        #ponteiro_flash_3phase2 = ponteiro_flash_3phase[ponteiro_flash_3phase]
        #ponteiro_flash_3phase2 = np.ones(len(self.P), dtype = bool)
        ponteiro_flash_3phase2[(sp2>1).sum(axis=0,dtype=bool)] = False
        #ponteiro_flash_3phase[ponteiro_flash_3phase] = ponteiro_flash_3phase2

        sp3, Kvalue3 = self.Stability_3phase(PR, self.y, np.copy(ponteiro_flash_3phase2))
        sp3 = np.round(sp3, 8)
        #import pdb; pdb.set_trace()
        ponteiro_aux2 = ~ponteiro_flash_3phase2[ponteiro_flash_3phase2]
        #ponteiro_aux2 = np.zeros(len(self.P), dtype = bool)
        ponteiro_aux2[(sp3>1).sum(axis=0,dtype=bool)] = True
        ponteiro_flash_3phase2[ponteiro_flash_3phase2] = ponteiro_aux2
        #ponteiro_flash_3phase2 = ponteiro_aux2

        ponteiro_flash_3phase = ponteiro_flash_3phase + ponteiro_flash_3phase2
        #ponteiro_flash_3phase[ponteiro_flash_3phase] = ponteiro_flash_3phase2
        import pdb; pdb.set_trace()



        '''
        ponteiro_flash_3phase2 = ~ponteiro_flash_3phase[ponteiro_flash_3phase]
        ponteiro_flash_3phase2[(sp3>1).sum(axis=0,dtype=bool)] = True

        #ponteiro_flash_3phase[ponteiro_flash_3phase] = ponteiro_flash_3phase2
        ponteiro_flash_3phase[ponteiro_flash_3phase2] = ponteiro_flash_3phase2
        '''

        'here!'
        self.K_A = self.K[:, ponteiro_flash_3phase]
        self.K_V = self.Kwilson[:, ponteiro_flash_3phase]

        self.molar_properties_3phase(PR, z, ponteiro_flash_3phase)
        self.get_other_properties_3phases(PR, Mw)
        import pdb; pdb.set_trace()
        ##########################
        if self.L != 1 and self.V != 1:
            sp2, Kvalue2 = self.Stability_3phase(PR, self.x, np.copy(ponteiro_flash))
            sp2 = np.round(sp2, 8)
            if all(sp2 <= 1):
                sp2, Kvalue2 = self.Stability_3phase(PR, self.y, np.copy(ponteiro_flash))

            sp2 = np.round(sp2, 8)
            print(f'sp2: {sp2}')
            #import pdb; pdb.set_trace()
            #print(f'K : {Kvalue2}')
            #sp2[0] = 5
            if any(sp2>1):
                index_sp2max = np.argmax(sp2)
                self.K_A = self.K.copy()
                #self.K_V = Kvalue2[1].copy()
                self.K_V = self.Kwilson.copy()
                self.molar_properties_3phase(PR, z, np.ones_like(ponteiro_flash, dtype=bool))
                self.get_other_properties_3phases(PR, Mw)
        ############################

        enthalpy_oil = PR.enthalpy_calculation(self, self.x, self.P, self.ph_L) # Teste de calculo da entalpia
        enthalpy_vapour = PR.enthalpy_calculation(self, self.y, self.P, self.ph_V) # Teste de calculo da entalpia

        import pdb; pdb.set_trace()
        #PR.get_all_derivatives(self, self.x, self.y, self.L, self.V, self.P)


    def get_other_properties(self, PR, Mw):
        self.Mw_L, self.ksi_L, self.rho_L = self.other_properties(PR, self.x, Mw, self.ph_L)
        self.Mw_V, self.ksi_V, self.rho_V = self.other_properties(PR, self.y, Mw, self.ph_V)

    def get_other_properties_3phases(self, PR, Mw):
        self.Mw_L, self.ksi_L, self.rho_L = self.other_properties(PR, self.x, Mw, self.ph_L)
        self.Mw_V, self.ksi_V, self.rho_V = self.other_properties(PR, self.y, Mw, self.ph_V)
        self.Mw_A, self.ksi_A, self.rho_A = self.other_properties(PR, self.a, Mw, self.ph_L)

    def equilibrium_ratio_Wilson(self):

        self.K = np.exp(5.37 * (1 + self.w) * (1 - 1 / (self.T / self.Tc)), dtype=np.double)[:,np.newaxis] / \
                (self.P / self.Pc[:,np.newaxis])
        self.Kwilson = self.K.copy()

    def equilibrium_ratio_aqueous(self, z):

        self.Kw = np.zeros_like(z)
        self.Kw[0] = 0.999 / z[0][0]
        self.Kw[1:] = 0.001 / (len(z) - 1) / z[1:,[0]]

        #self.Kw = 0.001 / (len(z) - 1) / z
        #self.Kw[3] = 0.999 / z[-1]


    def equilibrium_ratio_2flash(self, K_2flash):
        self.K = K_2flash.copy()

    def pressure_water_saturation(self):
        ' Correlation of Wagner and Saul 1987 '
        a1 = -7.85823
        a2 = 1.83991
        a3 = -11.7811
        a4 = 22.6705
        a5 = -15.9393
        a6 = 1.77516

        Tc = 647.14 # Kelvin
        Pc = 22.064e6 # Pascal
        tal = 1 - (self.T/Tc)

        ln_P_Pc = (Tc/self.T)*(a1*tal + a2*(tal**1.5) + a3*(tal**3) + a4*(tal**3.5) + a5*(tal**4) + a6*(tal**7.5))
        self.Pw_sat = np.exp(ln_P_Pc) * Pc

    """-------------Below starts stability test calculation -----------------"""
    """ delete this function """
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













    """-------------Below starts phase stability test calculation -----------------"""
    ''' delete this function '''
    def Stability_2phase(self, PR, z, ponteiro_stab_check):
        ''' In the lnphi function: 0 stands for vapor phase and 1 for liquid '''
        t0 = time.time()
    #*****************************Test one**********************************#
        stationary_points = np.empty((5, len(ponteiro_stab_check[ponteiro_stab_check])))
        # Trial phase is liquid
        Y = np.empty(z.shape)
        lnphiz = np.empty(z.shape)
        ponteiro = np.copy(ponteiro_stab_check)
        Y[:,ponteiro] = z[:,ponteiro] / self.K[:,ponteiro]
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphiz[:,ponteiro] = PR.lnphi(self, z[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
        K_value1 = 1 / self.K[:,ponteiro]

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = PR.lnphi(self, y[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(z[:,ponteiro]) + lnphiz[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        #lny = np.log(y)
        #TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphiz[:] - np.log(z)))
        #print(f'TPD do teste 1: {TPD}')
        #import pdb; pdb.set_trace()
        stationary_point1 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[0] = stationary_point1

    #*****************************Test two**********************************#
        # Trial phase is vapour
        ponteiro = np.copy(ponteiro_stab_check)

        Y[:,ponteiro] = self.K[:,ponteiro] * z[:,ponteiro]
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphiz[:,ponteiro] = PR.lnphi(self, z[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
        K_value2 = self.K[:,ponteiro]

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = PR.lnphi(self, y[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(z[:,ponteiro]) + lnphiz[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        lny = np.log(y)
        TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphiz[:] - np.log(z)))
        print(f'TPD do teste 2: {TPD}')

        stationary_point2 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[1] = stationary_point2


    #*****************************Test three**********************************#
        # Trial phase is liquid
        ponteiro = np.copy(ponteiro_stab_check)

        Y[:,ponteiro] = z[:,ponteiro] / ((self.K[:,ponteiro])**(1/3))
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphiz[:,ponteiro] = PR.lnphi(self, z[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
        K_value3 = 1 / ((self.K[:,ponteiro])**(1/3))

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = PR.lnphi(self, y[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(z[:,ponteiro]) + lnphiz[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        lny = np.log(y)
        TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphiz[:] - np.log(z)))
        print(f'TPD do teste 3: {TPD}')

        stationary_point3 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[2] = stationary_point3



    #*****************************Test four**********************************#
        # Trial phase is vapour
        ponteiro = np.copy(ponteiro_stab_check)

        Y[:,ponteiro] = z[:,ponteiro] * ((self.K[:,ponteiro])**(1/3))
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphiz[:,ponteiro] = PR.lnphi(self, z[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
        K_value4 = (self.K[:,ponteiro])**(1/3)

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = PR.lnphi(self, y[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(z[:,ponteiro]) + lnphiz[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        lny = np.log(y)
        TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphiz[:] - np.log(z)))
        print(f'TPD do teste 4: {TPD}')

        stationary_point4 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[3] = stationary_point4


    #*****************************Test five**********************************#
        # Trial phase is aqueous
        ponteiro = np.copy(ponteiro_stab_check)
        Y[:,ponteiro] = z[:,ponteiro] * self.Kw[:,ponteiro]
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]

        lnphiz[:,ponteiro] = PR.lnphi(self, z[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro]) # both ph_L ?
        K_value5 = self.Kw[:,ponteiro]

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = PR.lnphi(self, y[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro]) # both ph_L ?
            Y[:,ponteiro] = np.exp(np.log(z[:,ponteiro]) + lnphiz[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux


        lny = np.log(y)
        TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphiz[:] - np.log(z)))
        print(f'TPD do teste 5: {TPD}')

        stationary_point5 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[4] = stationary_point5

        t1 = time.time()
        print('stability test time:', t1-t0)

        #return stationary_point1, stationary_point2, stationary_point3, stationary_point4, stationary_point5
        return stationary_points, [K_value1, K_value2, K_value3, K_value4, K_value5]








    """-------------Below starts 3 phase stability test calculation -----------------"""

    def Stability_3phase(self, PR, x, ponteiro_stab_check):
        'Testing stability in system'
    #*****************************Test one**********************************#
        t0 = time.time()

        stationary_points = np.empty((5, len(ponteiro_stab_check[ponteiro_stab_check])))
        # Trial phase is liquid
        #import pdb; pdb.set_trace()
        Y = np.empty(x.shape)
        lnphix = np.empty(x.shape)
        ponteiro = np.copy(ponteiro_stab_check)
        Y[:,ponteiro] = x[:,ponteiro] / self.Kwilson[:,ponteiro]
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphix[:,ponteiro] = PR.lnphi(self, x[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
        K_value1 = 1 / self.Kwilson[:,ponteiro]

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = PR.lnphi(self, y[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(x[:,ponteiro]) + lnphix[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        #lny = np.log(y[:, ponteiro_stab_check])
        #TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphix[:] - np.log(x)))
        #print(f'TPD do teste 1: {TPD}')

        stationary_point1 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[0] = stationary_point1

    #*****************************Test two**********************************#
        # Trial phase is vapour
        ponteiro = np.copy(ponteiro_stab_check)

        Y[:,ponteiro] = self.Kwilson[:,ponteiro] * x[:,ponteiro]
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphix[:,ponteiro] = PR.lnphi(self, x[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
        K_value2 = self.Kwilson[:,ponteiro]

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = PR.lnphi(self, y[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(x[:,ponteiro]) + lnphix[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux


        lny = np.log(y)

        #TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphix[:] - np.log(x)))
        #print(f'TPD do teste 2: {TPD}')

        stationary_point2 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[1] = stationary_point2


    #*****************************Test three**********************************#
        # Trial phase is liquid
        ponteiro = np.copy(ponteiro_stab_check)

        Y[:,ponteiro] = x[:,ponteiro] / ((self.Kwilson[:,ponteiro])**(1/3))
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphix[:,ponteiro] = PR.lnphi(self, x[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
        K_value3 = 1 / ((self.Kwilson[:,ponteiro])**(1/3))

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = PR.lnphi(self, y[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(x[:,ponteiro]) + lnphix[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        lny = np.log(y)

        #TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphix[:] - np.log(x)))
        #print(f'TPD do teste 3: {TPD}')

        stationary_point3 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[2] = stationary_point3



    #*****************************Test four**********************************#
        # Trial phase is vapour
        ponteiro = np.copy(ponteiro_stab_check)

        Y[:,ponteiro] = x[:,ponteiro] * ((self.Kwilson[:,ponteiro])**(1/3))
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphix[:,ponteiro] = PR.lnphi(self, x[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
        K_value4 = (self.Kwilson[:,ponteiro])**(1/3)

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = PR.lnphi(self, y[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(x[:,ponteiro]) + lnphix[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        lny = np.log(y)

        #TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphix[:] - np.log(x)))
        #print(f'TPD do teste 4: {TPD}')

        stationary_point4 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[3] = stationary_point4


    #*****************************Test five**********************************#
        # Trial phase is aqueous
        ponteiro = np.copy(ponteiro_stab_check)
        Y[:,ponteiro] = x[:,ponteiro] * self.Kw[:,ponteiro]
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]

        lnphix[:,ponteiro] = PR.lnphi(self, x[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro]) # both ph_L ?
        K_value5 = self.Kw[:,ponteiro]

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = PR.lnphi(self, y[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro]) # both ph_L ?
            Y[:,ponteiro] = np.exp(np.log(x[:,ponteiro]) + lnphix[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        lny = np.log(y)

        #TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphix[:] - np.log(x)))
        #print(f'TPD do teste 5: {TPD}')

        stationary_point5 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[4] = stationary_point5

        #print(f'pontos estacionarios: {stationary_points}')

        t1 = time.time()
        print('stability time:', t1-t0)

        #return stationary_point1, stationary_point2, stationary_point3, stationary_point4, stationary_point5
        return stationary_points, [K_value1, K_value2, K_value3, K_value4, K_value5]












    """-------------Below starts biphasic flash calculations-----------------"""
    def molar_properties(self, PR, z, ponteiro):
        ponteiro = self.molar_properties_Yinghui(PR, z, ponteiro)
        #ponteiro = self.molar_properties_Whitson(PR, z, ponteiro)
        return ponteiro

    def molar_properties_3phase(self, PR, z, ponteiro):
        #ponteiro = self.molar_properties_Whitson_3phase(PR, z, ponteiro)
        ponteiro = self.molar_properties_Lapene_3phase(PR, z, ponteiro)
        return ponteiro

    def deltaG_molar_vectorized(self, PR, l, P, ph):
        lnphi = np.empty([2, len(self.w), len(ph)])

        try:
            lnphi[0,:] = PR.lnphi(self, l, P, 1 - ph)
        except: print('erro')#import pdb; pdb.set_trace()


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

    def solve_objective_function_Yinghui(self, z1, zi, K1, KNc, Ki, K, x, x1, i):
        x1_min = z1 * (1 - KNc) / (K1 - KNc)
        x1_max = (1 - KNc) / (K1 - KNc)

        vols_zi_neg = np.zeros(len(K1), dtype = bool)
        vols_zi_neg[np.sum(zi < 0, axis = 0, dtype=bool)] = True
        KNc_z_neg = KNc[vols_zi_neg].ravel()
        K1_z_neg = K1[vols_zi_neg]
        z1_z_neg = z1[vols_zi_neg]
        zi_z_neg = zi[:,vols_zi_neg]
        Ki_z_neg = Ki[:,vols_zi_neg]
        vols_zi_neg_num = np.sum(vols_zi_neg*1) + 1 - np.sign(np.sum(vols_zi_neg*1))
        Ki_z_neg_K_big1 = Ki_z_neg[Ki_z_neg > 1]#.reshape(int(len(Ki_z_neg[Ki_z_neg>1])/vols_zi_neg_num), vols_zi_neg_num)
        theta = np.ones(zi[:,vols_zi_neg].shape)

        theta[Ki_z_neg > 1] = ((1 - KNc_z_neg) / (Ki_z_neg_K_big1 - KNc_z_neg))#.ravel()
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

        x1_max_aux[cond_aux] = np.min(aux_eq_cond, axis = 0, initial=0)
        x1_max[vols_zi_neg] = x1_max_aux

        x1_min_aux = np.copy(x1_min[vols_zi_neg])
        x1_min_aux[~cond_aux] = np.max(aux_eq[~cond_aux], axis = 0, initial=0)

        x1_min_aux[x1_min_aux < 0] = 0
        x1_min[vols_zi_neg] = x1_min_aux

        if any(x1_min > x1_max):
            raise ValueError('There is no physical root')

        x1 = ((x1_min + x1_max) / 2)

        x1_new = np.copy(x1)
        ponteiro = np.ones(len(x1), dtype = bool)

        while any(ponteiro):
            x1[ponteiro] = np.copy(x1_new[ponteiro])
            f = 1 + ((K1[ponteiro] - KNc[ponteiro]) / (KNc[ponteiro] - 1)) * x1[ponteiro] + np.sum(((Ki[:,ponteiro] - KNc[ponteiro][np.newaxis,:]) /
                (KNc[ponteiro][np.newaxis,:] - 1)) * zi[:,ponteiro] * (K1[ponteiro][np.newaxis,:] - 1) * x1[ponteiro][np.newaxis,:]
                / ((Ki[:,ponteiro] - 1) * z1[ponteiro][np.newaxis,:] + (K1[ponteiro][np.newaxis,:] - Ki[:,ponteiro]) *
                x1[ponteiro][np.newaxis,:]), axis = 0)
            df = ((K1[ponteiro] - KNc[ponteiro]) / (KNc[ponteiro] - 1)) + np.sum(((Ki[:,ponteiro] - KNc[ponteiro][np.newaxis,:]) /
                (KNc[ponteiro][np.newaxis,:] - 1)) * zi[:,ponteiro] * z1[ponteiro][np.newaxis,:] * (K1[ponteiro][np.newaxis,:] - 1) *
                (Ki[:,ponteiro] - 1) / ((Ki[:,ponteiro] - 1) * z1[ponteiro][np.newaxis,:] + (K1[ponteiro][np.newaxis,:] - Ki[:,ponteiro]) *
                x1[ponteiro][np.newaxis,:]) ** 2, axis = 0)
            x1_new[ponteiro] = x1[ponteiro] - f/df #Newton-Raphson iterative method
            x1_aux = x1_new[ponteiro]
            x1_aux[x1_aux > x1_max] = (x1_min[x1_aux > x1_max] + x1_max[x1_aux > x1_max])/2
            x1_aux[x1_aux < x1_min] = (x1_min[x1_aux < x1_min] + x1_max[x1_aux < x1_min])/2
            x1_new[ponteiro] = x1_aux
            ponteiro_aux = ponteiro[ponteiro] #o que muda de tamanho
            ponteiro_aux[abs(f) < 1e-10] = False
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


    def Yinghui_method(self, z, ponteiro, x1, i):

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
                                            K[:,~(z1 == 0)], x[:,~(z1 == 0)], x1[~(z1 == 0)], i)


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
        t0 = time.time()

        i=0
        x1 = np.zeros_like(z[0])
        while any(ponteiro):
            i+=1
            self.Yinghui_method(z[:,ponteiro], ponteiro, x1[ponteiro], i)
            lnphil = PR.lnphi(self, self.x[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            lnphiv = PR.lnphi(self, self.y[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            self.fl = np.exp(lnphil) * (self.x[:,ponteiro] * self.P[ponteiro][np.newaxis,:])
            self.fv = np.exp(lnphiv) * (self.y[:,ponteiro] * self.P[ponteiro][np.newaxis,:])
            razao[:,ponteiro] = np.divide(self.fl, self.fv, out = razao[:,ponteiro] / razao[:,ponteiro] * (1 + 1e-10),
                              where = self.fv != 0)
            #x1 = self.x[self.K==np.max(self.K,axis=0)]
            self.K[:,ponteiro] = razao[:,ponteiro] * self.K[:,ponteiro]
            stop_criteria = np.max(abs(razao[:,ponteiro] - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        print('Yinghui it:', i)

        V = (z[:,ponteiro_save][self.x[:,ponteiro_save] != 0] - self.x[:,ponteiro_save][self.x[:,ponteiro_save] != 0]) / \
                          (self.y[:,ponteiro_save][self.x[:,ponteiro_save] != 0] - self.x[:,ponteiro_save][self.x[:,ponteiro_save] != 0])
        #vols_V = np.sum(self.x[:,ponteiro_save] == 0, dtype=bool, axis = 0)
        vv = np.argwhere((self.x[:,ponteiro_save]!=0) == True)
        vols_V, ind = np.unique(vv[:,1],return_index = True)
        self.V[ponteiro_save] = V[ind]
        self.L[ponteiro_save] = 1. - self.V[ponteiro_save]

        t1 = time.time()
        print('Yinghui time:', t1-t0)
        return ponteiro_save

    def solve_objective_function_Whitson_for_L(self, z, L, Lmax, Lmin, ponteiro):

        #Lmax = 1 / (1 - np.min(self.K, axis=0))
        #Lmin = ((np.max(self.K, axis=0) - np.min(self.K, axis=0))*z[self.K == np.max(self.K, axis = 0)] - (1 - np.min(self.K, axis=0))) / ((1 - np.min(self.K, axis=0))*(np.max(self.K, axis=0) - 1))
        #import pdb; pdb.set_trace()

        ponteiro_save = np.copy(ponteiro)
        while any(ponteiro):
             Lold = L[ponteiro]
             f = np.sum((1 - self.K[:,ponteiro]) * z[:,ponteiro] / (L[ponteiro][np.newaxis,:]
                        + (1 - L[ponteiro][np.newaxis,:]) * self.K[:,ponteiro]), axis = 0)
             df = - np.sum((1 - self.K[:,ponteiro]) ** 2 * z[:,ponteiro] / (L[ponteiro][np.newaxis,:] +
                    (1 - L[ponteiro][np.newaxis,:]) * self.K[:,ponteiro]) ** 2, axis = 0)
             L[ponteiro] = L[ponteiro] - f / df #Newton-Raphson iterative method
             L_aux = L[ponteiro]
             L_aux[L_aux > Lmax[ponteiro]] = (Lmax[ponteiro] + Lold)[L_aux > Lmax[ponteiro]] / 2
             L_aux[L_aux < Lmin[ponteiro]] = (Lmin[ponteiro] + Lold)[L_aux < Lmin[ponteiro]] / 2
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


        Vmax = 1 / (1 - np.min(self.K, axis=0))
        Vmin = ((np.max(self.K, axis=0) - np.min(self.K, axis=0))*z[self.K == np.max(self.K, axis = 0)] - (1 - np.min(self.K, axis=0))) / ((1 - np.min(self.K, axis=0))*(np.max(self.K, axis=0) - 1))
        #import pdb; pdb.set_trace()

        ponteiro_save = np.copy(ponteiro)
        i = 0
        while any(ponteiro):
            Vold = V[ponteiro]
            f = np.sum((self.K[:,ponteiro] - 1) * z[:,ponteiro] / (1 + V[ponteiro][np.newaxis,:] * (self.K[:,ponteiro] - 1)), axis = 0)
            df = - np.sum((self.K[:,ponteiro] - 1) ** 2 * z[:,ponteiro] / (1 + V[ponteiro][np.newaxis,:] * (self.K[:,ponteiro] - 1)) ** 2, axis = 0)
            V[ponteiro] = V[ponteiro] - f / df #Newton-Raphson iterative method
            V_aux = self.V[ponteiro]
            V_aux[V_aux > Vmax[ponteiro]] = 0.5 * (Vmax[ponteiro] + Vold)[V_aux > Vmax[ponteiro]] #(Vmax + Vold)/2
            V_aux[V_aux < Vmin[ponteiro]] = 0.5 * (Vmin[ponteiro] + Vold)[V_aux < Vmin[ponteiro]]#(Vmax + Vold)/2
            #V_aux[V_aux > Vmax[ponteiro]] = Vmax[ponteiro][V_aux > Vmax[ponteiro]]
            #V_aux[V_aux < Vmin[ponteiro]] = Vmin[ponteiro][V_aux < Vmin[ponteiro]]
            self.V[ponteiro] = V_aux
            stop_criteria = abs(V[ponteiro] / Vold - 1)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux
            i+=1
            if i>=100:
                print('maxit')
                return False

            #import pdb; pdb.set_trace()
        #self.V[ponteiro_save] = V[ponteiro_save] + 1e-15 #manipulação
        self.x[:,ponteiro_save] = z[:,ponteiro_save] / (1 + self.V[ponteiro_save][np.newaxis,:] * (self.K[:,ponteiro_save] - 1))
        self.y[:,ponteiro_save] = self.K[:,ponteiro_save] * self.x[:,ponteiro_save]
        self.L[ponteiro_save] = 1. - self.V[ponteiro_save]

        return True


    def molar_properties_Whitson(self, PR, z, ponteiro):
        Lmax = np.max(self.K, axis = 0)/(np.max(self.K, axis = 0) - 1)
        Lmin = np.min(self.K, axis = 0)/(np.min(self.K, axis = 0) - 1)
        Vmax = 1. - Lmin
        Vmin = 1. - Lmax
        #Vmin = ((K1-KNc)z[self.K==K1]-(1-KNc))/((1-KNc)(K1-1))
        #proposed by Li et al for Whitson method
        self.V[ponteiro] = (Vmin[ponteiro] + Vmax[ponteiro]) * 0.5
        self.L = 1 - self.V
        t0 = time.time()

        razao = np.ones(z.shape)/2
        ponteiro_save = np.copy(ponteiro)
        while any(ponteiro):
            #try_V = self.solve_objective_function_Whitson_for_V(z, self.V, Vmax, Vmin, np.copy(ponteiro))
            #if not try_V:
                #print('Entrou no L')
            #self.solve_objective_function_Whitson_for_V(z, self.V, Vmax, Vmin, np.copy(ponteiro))
            self.solve_objective_function_Whitson_for_L(z, self.L, Lmax, Lmin, np.copy(ponteiro))
            #else: self.L[ponteiro] = 1. - self.V[ponteiro]
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
            ponteiro[((self.V)<0) + ((self.V)>1)] = False

        t1 = time.time()
        print('Whitson-Michelsen time:', t1-t0)

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

        #import pdb; pdb.set_trace()
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








    ' Try of triphasic flash '
    def molar_properties_Whitson_3phase(self, PR, z, ponteiro):
        Amax = 1 / (1 - np.min(self.K_A, axis = 0))
        Amin = 1 / (1 - np.max(self.K_A, axis = 0))
        Vmax = 1 / (1 - np.min(self.K_V, axis = 0))
        Vmin = 1 / (1 - np.max(self.K_V, axis = 0))

        Amax = np.max(self.K_A, axis = 0)/(np.max(self.K_A, axis = 0) - 1)
        Amin = np.min(self.K_A, axis = 0)/(np.min(self.K_A, axis = 0) - 1)

        #Amax = 1. - Vmin
        #Amin = 1. - Vmax

        #Amin = np.array([0.0])
        #Vmax = np.array([1.0])
        #Vmin = np.array([0.0])

        self.V[ponteiro] = (Vmin[ponteiro] + Vmax[ponteiro]) * 0.5
        self.A[ponteiro] = (Amin[ponteiro] + Amax[ponteiro]) * 0.5

        t0 = time.time()

        razao_vl = np.ones(z.shape)/2
        razao_al = np.ones(z.shape)/2
        ponteiro_save = np.copy(ponteiro)
        while any(ponteiro):
            self.solve_objective_function_Whitson_for_V_and_A(z, self.V, self.A, Vmax, Vmin, Amax, Amin, np.copy(ponteiro))

            lnphil = self.lnphi_based_on_deltaG(PR, self.x[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            lnphiv = self.lnphi_based_on_deltaG(PR, self.y[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            lnphia = self.lnphi_based_on_deltaG(PR, self.a[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])

            self.fv = np.exp(lnphiv) * (self.y[:,ponteiro] * self.P[ponteiro][np.newaxis,:])
            self.fl = np.exp(lnphil) * (self.x[:,ponteiro] * self.P[ponteiro][np.newaxis,:])
            self.fa = np.exp(lnphia) * (self.a[:,ponteiro] * self.P[ponteiro][np.newaxis,:])

            razao_vl[:,ponteiro] = np.divide(self.fl, self.fv, out = razao_vl[:,ponteiro] / razao_vl[:,ponteiro] * (1 + 1e-10),
                              where = self.fv != 0)
            razao_al[:,ponteiro] = np.divide(self.fl, self.fa, out = razao_al[:,ponteiro] / razao_al[:,ponteiro] * (1 + 1e-10),
                              where = self.fa != 0)

            self.K_V[:,ponteiro] = razao_vl[:,ponteiro] * self.K_V[:,ponteiro]
            self.K_A[:,ponteiro] = razao_al[:,ponteiro] * self.K_A[:,ponteiro]

            stop_criteria_vl = np.max(abs((self.fv/(self.fl + 1e-15)) - 1), axis = 0)
            stop_criteria_al = np.max(abs((self.fa/(self.fl + 1e-15)) - 1), axis = 0)
            stop_criteria = max(stop_criteria_vl, stop_criteria_al)

            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux
            #ponteiro[((self.V)<0) + ((self.V)>1)] = False

        t1 = time.time()
        print('Whitson-Michelsen time for 3phase flash:', t1-t0)

        return ponteiro_save

    def solve_objective_function_Whitson_for_V_and_A(self, z, V, A, Vmax, Vmin, Amax, Amin, ponteiro):

        ponteiro_save = np.copy(ponteiro)
        i = 0
        V_and_A = np.array([V, A])
        while any(ponteiro):
            Vold = V[ponteiro].copy()
            Aold = A[ponteiro].copy()
            #V_and_A = np.array([Vold, Aold])
            #print(V_and_A)

            f = np.zeros(2)

            f[0] = np.sum((self.K_V[:,ponteiro] - 1) * z[:,ponteiro] / (1 + \
                    V[ponteiro][np.newaxis,:] * (self.K_V[:,ponteiro] - 1) + \
                    A[ponteiro][np.newaxis,:] * (self.K_A[:,ponteiro] - 1)), axis = 0)

            f[1] = np.sum((self.K_A[:,ponteiro] - 1) * z[:,ponteiro] / (1 + \
                    V[ponteiro][np.newaxis,:] * (self.K_V[:,ponteiro] - 1) + \
                    A[ponteiro][np.newaxis,:] * (self.K_A[:,ponteiro] - 1)), axis = 0)

            df_VdV = - np.sum(((self.K_V[:,ponteiro] - 1) ** 2) * z[:,ponteiro] / ((1 + \
                    V[ponteiro][np.newaxis,:] * (self.K_V[:,ponteiro] - 1) + \
                    A[ponteiro][np.newaxis,:] * (self.K_A[:,ponteiro] - 1)) ** 2), axis = 0)

            df_VdA = - np.sum((self.K_V[:,ponteiro] - 1) * (self.K_A[:,ponteiro] - 1) * z[:,ponteiro] / ((1 + \
                    V[ponteiro][np.newaxis,:] * (self.K_V[:,ponteiro] - 1) + \
                    A[ponteiro][np.newaxis,:] * (self.K_A[:,ponteiro] - 1)) ** 2), axis = 0)

            df_AdV = - np.sum((self.K_V[:,ponteiro] - 1) * (self.K_A[:,ponteiro] - 1) * z[:,ponteiro] / ((1 + \
                    V[ponteiro][np.newaxis,:] * (self.K_V[:,ponteiro] - 1) + \
                    A[ponteiro][np.newaxis,:] * (self.K_A[:,ponteiro] - 1)) ** 2), axis = 0)

            df_AdA = - np.sum(((self.K_A[:,ponteiro] - 1) ** 2) * z[:,ponteiro] / ((1 + \
                    V[ponteiro][np.newaxis,:] * (self.K_V[:,ponteiro] - 1) + \
                    A[ponteiro][np.newaxis,:] * (self.K_A[:,ponteiro] - 1)) ** 2), axis = 0)

            Jacobian = np.empty([2,2])
            Jacobian[0][0] = df_VdV
            Jacobian[0][1] = df_VdA
            Jacobian[1][0] = df_AdV
            Jacobian[1][1] = df_AdA # Jacobian matrix

            try:
                Jacobian_inv = np.linalg.inv(Jacobian)
            except: import pdb; pdb.set_trace()

            V_and_A = V_and_A - np.reshape(f.dot(Jacobian_inv), (2,1))  #Newton-Raphson iterative method
            V = V_and_A[0].copy()
            A = V_and_A[1].copy()

            V_aux = V[ponteiro]
            V_aux[V_aux > Vmax[ponteiro]] = 0.5 * (Vmax[ponteiro] + Vold)[V_aux > Vmax[ponteiro]] #(Vmax + Vold)/2
            V_aux[V_aux < Vmin[ponteiro]] = 0.5 * (Vmin[ponteiro] + Vold)[V_aux < Vmin[ponteiro]]#(Vmin + Vold)/2
            V[ponteiro] = V_aux

            A_aux = A[ponteiro]
            A_aux[A_aux > Amax[ponteiro]] = 0.5 * (Amax[ponteiro] + Aold)[A_aux > Amax[ponteiro]] #(Amax + Aold)/2
            A_aux[A_aux < Amin[ponteiro]] = 0.5 * (Amin[ponteiro] + Aold)[A_aux < Amin[ponteiro]]#(Amin + Aold)/2
            A[ponteiro] = A_aux

            stop_criteria = max(abs((V[ponteiro] / Vold) - 1), abs((A[ponteiro] / Aold) - 1))
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux
            i+=1
            if i>=3000:
                print('maxit in triphasic flash')
                import pdb; pdb.set_trace()

        self.V[ponteiro_save] = V[ponteiro_save]
        self.A[ponteiro_save] = A[ponteiro_save]
        self.L[ponteiro_save] = (1 - self.V[ponteiro_save] - self.A[ponteiro_save])

        self.x[:,ponteiro_save] = z[:,ponteiro_save] / (1 + self.V[ponteiro_save][np.newaxis,:] * (self.K_V[:,ponteiro_save] - 1) + \
                                    self.A[ponteiro_save][np.newaxis,:] * (self.K_A[:,ponteiro_save] - 1))
        self.y[:,ponteiro_save] = self.K_V[:,ponteiro_save] * self.x[:,ponteiro_save]
        self.a[:,ponteiro_save] = self.K_A[:,ponteiro_save] * self.x[:,ponteiro_save]

        return True




    def molar_properties_Lapene_3phase(self, PR, z, ponteiro):

        V = self.V[ponteiro]
        x = self.x[:, ponteiro]
        y = self.y[:, ponteiro]

        self.K_V[0] = (self.Pc[0]/self.P[ponteiro])*(self.T/self.Tc[0]) # Adicionar ponteiro no T para o termico

        y[0] = self.Pw_sat / self.P[ponteiro]
        self.K2w = y[0].copy()
        x[0] = y[0] / self.K_V[0]
        self.a[:, ponteiro] = 0
        self.a[0, ponteiro] = 1

        t0 = time.time()

        razao_vl = np.ones(x.shape)/2
        razao_av = np.ones(self.P[ponteiro].shape)/2
        ponteiro_save = np.copy(ponteiro)

        while any(ponteiro):
            y[0] = self.K2w.copy()
            x[0] = y[0] / self.K_V[0]

            K_V_max = np.amax(self.K_V[1:,:], axis=0)
            K_V_min = np.amin(self.K_V[1:,:], axis=0)
            Kw_ast = (1 - y[0]) / (1 - x[0])
            Kwz = (1 - z[0, ponteiro]) / (1 - x[0])

            Vmax = Kwz / (Kw_ast - K_V_min)
            Vmin = Kwz / (Kw_ast - K_V_max)

            ponteiro[ponteiro] = ~((K_V_max < Kw_ast) + (K_V_min > Kw_ast))
            '''
            if (K_V_max < Kw_ast or K_V_min > Kw_ast):
                print('Solve for V a classical RR equation')
                # colocar False no ponteiro
                import pdb; pdb.set_trace()
            '''

            V_dash = (z[0, ponteiro] - x[0]) / (y[0] - x[0])
            Objetive_function_V_dash = np.sum((self.K_V[1:] - Kw_ast) * z[1:,ponteiro] / (Kwz + \
                                        V_dash * (self.K_V[1:] - Kw_ast)), axis = 0)

            # Duvida se o calculo de cima esta certo mesmo

            ' determine position of V with respect V_dash '
            V_bigger_then_Vdash = np.ones(len(V), dtype = bool)
            ponteiro_aux = np.ones(len(V), dtype = bool)

            ponteiro_aux[ponteiro_aux] = (V_dash > Vmin) & (V_dash < Vmax)
            V_bigger_then_Vdash[ponteiro_aux] = Objetive_function_V_dash[ponteiro_aux] > 0
            '''
            if (V_dash > Vmin and V_dash < Vmax):
                if (Objetive_function_V_dash > 0):
                    V_bigger_then_Vdash = True
                else:
                    V_bigger_then_Vdash = False
            '''

            #ponteiro_aux[~ponteiro_aux] = (V_dash < Vmin)[~ponteiro_aux]
            ponteiro_aux = V_dash < Vmin
            V_bigger_then_Vdash[ponteiro_aux] = True
            '''
            elif (V_dash < Vmin):
                V_bigger_then_Vdash = True
            '''

            #ponteiro_aux[~ponteiro_aux] = V_dash > Vmax
            ponteiro_aux = V_dash > Vmax
            V_bigger_then_Vdash[ponteiro_aux] = False
            '''
            elif (V_dash > Vmax):
                V_bigger_then_Vdash = False
            '''

            ' determine presence or absence of water '
            cond_ext = ~V_bigger_then_Vdash
            cond_int = y[0] < x[0]

            cond1 = cond_ext & cond_int
            cond2 = ~cond_ext & ~cond_int

            aux = self.A[ponteiro]
            aux[cond1 + cond2] = 0.0
            self.A[ponteiro] = aux

            aux2 = ponteiro[ponteiro]
            aux2[cond1 + cond2] = False
            ponteiro[ponteiro] = aux2
            import pdb; pdb.set_trace()
            '''
            if not V_bigger_then_Vdash:
                if (y[0] < x[0]):
                    self.A = 0
                    print('Não tem água no sistema')
                    # ponteiro falso
                else:
                    print('Tem água no sistema')
            if V_bigger_then_Vdash:
                if (y[0] > x[0]):
                    self.A = 0
                    print('Não tem água no sistema')
                    # ponteiro falso
                else:
                    print('Tem água no sistema')
            '''

            """---------------- Stop Here ------------"""
            V, x, y = self.solve_objective_function_Lapene_for_V(z, x, y, V, Kw_ast, Kwz, np.copy(ponteiro))

            lnphil = self.lnphi_based_on_deltaG(PR, x[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            lnphiv = self.lnphi_based_on_deltaG(PR, y[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            lnphia = self.lnphi_based_on_deltaG(PR, self.a[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])


            self.fv = np.exp(lnphiv) * (y[:,ponteiro] * self.P[ponteiro][np.newaxis,:])
            self.fl = np.exp(lnphil) * (x[:,ponteiro] * self.P[ponteiro][np.newaxis,:])
            self.fa = np.exp(lnphia) * (self.a[:,ponteiro] * self.P[ponteiro][np.newaxis,:])


            razao_vl[:,ponteiro] = np.divide(self.fl, self.fv, out = razao_vl[:,ponteiro] / razao_vl[:,ponteiro] * (1 + 1e-10),
                              where = self.fv != 0)
            razao_av[ponteiro] = self.fa[0,ponteiro] / self.fv[0,ponteiro] * (1 + 1e-10)


            self.K_V[:,ponteiro] = razao_vl[:,ponteiro] * self.K_V[:,ponteiro]
            self.K2w = self.K2w * razao_av

            #stop_criteria = np.max(abs((self.fl/(self.fv + 1e-15)) - 1), axis = 0)
            stop_criteria = (np.sum(((self.fl/self.fv) - 1)**2))**0.5
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

            #ponteiro[((self.V)<0) + ((self.V)>1)] = False

        t1 = time.time()
        print('Lapene time for 3phase flash:', t1-t0)

        V_ast = (1 - z[0]) / (1 - y[0])
        if V > V_ast:
            V = V_ast
            print('V maior que o V asterisco')

        if V < 0:
            V = 0

        # Aqui atualizar tudo
        #self.V = V[ponteiro_save].copy()
        #self.x[:,ponteiro_save] = x.copy()
        #self.y[:,ponteiro_save] = y.copy()
        self.V = V
        self.x = x
        self.y = y
        self.A = (z[0] + self.V*(self.x[0] - self.y[0]) - self.x[0]) / (1 - self.x[0])
        self.L = 1 - self.A - self.V

        return ponteiro_save

    def solve_objective_function_Lapene_for_V(self, z, x, y, V, Kw_ast, Kwz, ponteiro):
        ' Leibovici–Neoschil window '
        Vmin = np.max(((self.K_V[self.K_V>Kw_ast] * z[self.K_V>Kw_ast] - Kwz) / (self.K_V[self.K_V>Kw_ast] - Kw_ast)))
        Vmax = np.min(((Kwz - z[self.K_V<Kw_ast]) / (Kw_ast - self.K_V[self.K_V<Kw_ast])), axis = 0)
        #problema com o ponteiro, Vmin e Vmax nao sao arrays
        V[ponteiro] = 0.5*(Vmax + Vmin) # Chute inicial

        ponteiro_save = np.copy(ponteiro)
        i = 0

        while any(ponteiro):
            Vold = V[ponteiro]

            G = np.sum((self.K_V[1:,ponteiro] - Kw_ast) * z[1:,ponteiro] / (Kwz + \
                        V * (self.K_V[1:,ponteiro] - Kw_ast)), axis = 0)

            dG = - np.sum(((self.K_V[1:,ponteiro] - Kw_ast)**2) * z[1:,ponteiro] / ((Kwz + \
                        V * (self.K_V[1:,ponteiro] - Kw_ast))**2), axis = 0)

            V[ponteiro] = V[ponteiro] - G / dG # Newton-Raphson iterative method

            V_aux = V[ponteiro]
            #V_aux[V_aux > Vmax[ponteiro]] = 0.5 * (Vmax[ponteiro] + Vold)[V_aux > Vmax[ponteiro]] #(Vmax + Vold)/2
            #V_aux[V_aux < Vmin[ponteiro]] = 0.5 * (Vmin[ponteiro] + Vold)[V_aux < Vmin[ponteiro]]#(Vmin + Vold)/2
            V_aux[V_aux > Vmax] = 0.5 * (Vmax + Vold)[V_aux > Vmax] #(Vmax + Vold)/2
            V_aux[V_aux < Vmin] = 0.5 * (Vmin + Vold)[V_aux < Vmin]#(Vmin + Vold)/2
            V[ponteiro] = V_aux

            stop_criteria = abs((V[ponteiro] / Vold) - 1)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux
            i+=1
            if i>=3000:
                print('maxit in triphasic flash')
                import pdb; pdb.set_trace()

        x[1:] = z[1:] / (1 + V*(self.K_V[1:] - 1 + (y[0] - x[0])/(1 - x[0])) + (x[0] - z[0])/(1 - x[0]))
        y[1:] = x[1:] * self.K_V[1:]
        y[0] = self.K2w.copy()
        x[0] = y[0] / self.K_V[0]

        #self.V[ponteiro_save] = V[ponteiro_save]
        #self.x[:,ponteiro_save] = x
        #self.y[:,ponteiro_save] = y
        return V, x, y

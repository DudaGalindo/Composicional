import numpy as np
from math import pi
from solver import CubicRoots
from scipy import linalg
import matplotlib.pyplot as plt

class PengRobinson:
    def __init__(self, kprop):

        self.coefficientsPR(kprop)

    def coefficientsPR(self, kprop):

        #l - any phase molar composition
        PR_kC7 = np.array([0.379642, 1.48503, 0.1644, 0.016667])
        PR_k = np.array([0.37464, 1.54226, 0.26992])

        k = (PR_kC7[0] + PR_kC7[1] * kprop.w - PR_kC7[2] * kprop.w ** 2 + \
            PR_kC7[3] * kprop.w ** 3) * (1*(kprop.w >= 0.5)) + (PR_k[0] + PR_k[1] * kprop.w - \
            PR_k[2] * kprop.w ** 2) * (1*(kprop.w < 0.5))  # FUnção com o fator acentrico - 2 possibilidades, varia com o w
        alpha = (1 + k * (1 - (kprop.T / kprop.Tc) ** (1 / 2))) ** 2
        aalpha_i = 0.45724 * (kprop.R * kprop.Tc) ** 2 / kprop.Pc * alpha
        self.b = 0.07780 * kprop.R * kprop.Tc / kprop.Pc
        aalpha_i_reshape = np.ones((kprop.Nc,kprop.Nc)) * aalpha_i[:,np.newaxis]
        self.aalpha_ik = np.sqrt(aalpha_i_reshape.T * aalpha_i[:,np.newaxis]) \
                        * (1 - kprop.Bin)

    def coefficients_cubic_EOS_all(self, kprop, l, P):
        self.bm = np.sum(l * self.b[:,np.newaxis], axis=0).ravel()
        l_reshape = np.ones((self.aalpha_ik).shape)[:,:,np.newaxis] * l[:,np.newaxis,:]
        #dd = l #tentar melhorar essa função
        #import pdb; pdb.set_trace()
        #self.aalpha = l @ l_reshape
        self.aalpha = (l_reshape * l[np.newaxis,:,:] * self.aalpha_ik[:,:,np.newaxis]).sum(axis=0).sum(axis=0)
        B = self.bm * P / (kprop.R* kprop.T)
        A = self.aalpha * P / (kprop.R* kprop.T) ** 2
        self.psi = (l_reshape * self.aalpha_ik[:,:,np.newaxis]).sum(axis = 0)
        return A, B

    def Z_all(B, A, ph):
        coef = np.array([np.ones(len(B)), -(1 - B), (A - 2*B - 3*B**2), -(A*B - B**2 - B**3)])

        Q = 2*coef[1]**3/27 - coef[1]*coef[2]/3 + coef[3]
        P = -coef[1]**2/3 + coef[2]
        delta = (Q/2)**2 + (P/3)**3
        R = (-P/3)**(3/2)
        theta = 1/3*np.arccos(-Q/(2*R))
        omega = (-1 + 1j*np.sqrt(3))/2

        X1 = np.cbrt(-Q/2 + delta**(1/2)) + np.cbrt(-Q/2 - delta**(1/2))
        X2 = omega*np.cbrt(-Q/2 + delta**(1/2)) + omega**2*np.cbrt(-Q/2 - delta**(1/2))
        X3 = omega**2*np.cbrt(-Q/2 + delta**(1/2)) + omega**3*np.cbrt(-Q/2 - delta**(1/2))
        X = np.array([X1,X2,X3])
        reais = np.isreal(X)
        r_pos = np.where(reais==True)
        Xreais = np.real(X[r_pos[:]])
        Z = Xreais -coef[1]/3
        Z_ans = np.min(Z, axis=0) * ph + np.max(Z, axis=0) * (1 - ph)
        return Z_ans
        #Z = np.linalg.lstsq(coef, np.zeros(len(A[:])))

    def coefficients_cubic_EOS(self, kprop, l, P):
        self.bm = sum(l * self.b)
        l_reshape = np.ones((self.aalpha_ik).shape) * l[:, np.newaxis]
        self.aalpha = (l_reshape.T * l[:,np.newaxis] * self.aalpha_ik).sum()
        B = self.bm * P / (kprop.R* kprop.T)
        A = self.aalpha * P / (kprop.R* kprop.T) ** 2
        self.psi = (l_reshape * self.aalpha_ik).sum(axis = 0)
        return A, B

    def Z(B, A, ph):
        # PR cubic EOS: Z**3 - (1-B)*Z**2 + (A-2*B-3*B**2)*Z-(A*B-B**2-B**3)
        coef = np.array([1, -(1 - B), (A - 2*B - 3*B**2), -(A*B - B**2 - B**3)])

        Z = np.roots(coef)
        root = np.isreal(Z) # return True for real roots
        #position where the real roots are - crated for organization only
        real_roots_position = np.where(root == True)
        Z_reais = np.real(Z[real_roots_position[:]]) #Saving the real roots
        Z_reais = Z_reais[Z_reais>0]
        Z_ans = min(Z_reais) * ph + max(Z_reais) * (1 - ph)

        ''' This last line, considers that the phase is composed by a pure
         component, so the EOS model can return more than one real root.
            If liquid, Zl = min(Z) and gas, Zv = max(Z).
            You can notice that, if there's only one real root,
        it works as well.'''
        return Z_ans

    def Z_vectorized(self, A, B, ph):
        coef = np.empty([4,len(B.ravel())])
        coef[0,:] = 1
        coef[1,:] = -(1 - B)
        coef[2,:] = (A - 2*B - 3*B**2)
        coef[3,:] = -(A*B - B**2 - B**3)
        Z = CubicRoots().run(coef)
        root = np.isreal(Z)
        n_reais = np.sum(root*1, axis = 1)

        aux_reais = (n_reais<3) & (n_reais>1)
        Z_save = Z
        if any(aux_reais): Z[~root[aux_reais]] = Z[root[aux_reais]][0]
        #Z[~root[n_reais>1]] = Z[root[n_reais > 1]][np.any(root[n_reais>1]==True, axis = 1)].ravel()
        Z[~root[n_reais==1]] = np.repeat(Z[root[n_reais == 1]], 2)
        aux_neg = np.zeros(Z.shape,dtype=bool)
        aux_neg[Z<0] = True

        try:
            Z[aux_neg] = Z[~aux_neg][0]
        except: import pdb; pdb.set_trace()
        #Z[aux_neg] = Z[~aux_neg][0]

        Zz = np.min(Z, axis = 1) * ph + np.max(Z, axis = 1) * (1 - ph)
        Z_ = np.real(Zz)
        #import pdb; pdb.set_trace()
        return Z_

    def lnphi(self, kprop, l, P, ph):
        #l - any phase molar composition
        A, B = self.coefficients_cubic_EOS_all(kprop,l, P)
        #Zfunc = np.vectorize(PengRobinson.Z)
        Z = self.Z_vectorized(A, B, ph)
        lnphi = self.lnphi_calculation(A, B, Z)
        return lnphi

    def lnphi_calculation(self, A, B, Z):
        #if Z==B: Z +=1e-20
        #Z = PengRobinson.Z_all(B, A, ph) #Not working here, but in the original code its working fine
        lnphi = self.b[:,np.newaxis] / self.bm * (Z[np.newaxis,:] - 1) - np.log(abs(Z[np.newaxis,:] - B[np.newaxis,:])) - A[np.newaxis,:] / (2 * (2 ** (1/2))
                * B[np.newaxis,:]) * (2 * self.psi / self.aalpha[np.newaxis,:] - self.b[:,np.newaxis] / self.bm) * np.log((Z[np.newaxis,:] + (1 +
                2 ** (1/2)) * B[np.newaxis,:]) / (Z[np.newaxis,:] + (1 - 2 ** (1/2)) * B[np.newaxis,:]))

        return lnphi

    def get_all_derivatives(self, kprop, x, y, L, V, P):

        #x = fprop.component_molar_fractions[0:kprop.Nc,0,:]
        #y = fprop.component_molar_fractions[0:kprop.Nc,1,:]
        #Nl = fprop.phase_mole_numbers[0,0,:]
        #Nv = fprop.phase_mole_numbers[0,1,:]
        #L[L>1] = 1
        #V[V<0] = 0
        Sw = 0.2
        Sg = (1. - Sw) * (V / kprop.ksi_V) / (V / kprop.ksi_V + L / kprop.ksi_L)
        So = 1 - Sg - Sw
        So = So/(So+Sg)
        Sg = Sg/(So+Sg)
        Vt = 11.32673864 #anything
        Nl = kprop.ksi_L * Vt * So
        Nv = kprop.ksi_V * Vt * Sg
        import pdb; pdb.set_trace()
        # REVER CONTAS
        dfildP, dfildnij, dZldP_parcial, dZldnij_parcial, Zl = \
                self.get_phase_derivatives(kprop, P, x, Nl, 1)
        dfivdP, dfivdnij, dZvdP_parcial, dZvdnij_parcial, Zv = \
                self.get_phase_derivatives(kprop, P, y, Nv, 0)

        dnldP, dnvdP, dnldNk, dnvdNk, dnildP, dnivdP, dnildNk, dnivdNk = self.dnij_dNk_dP(kprop, dfildP, dfivdP, dfildnij, dfivdnij, Nl, Nv)
        dZldP, dZvdP, dZldNk, dZvdNk = self.dZ_dP_dNk(dZldP_parcial, dZvdP_parcial,
                dZldnij_parcial, dZvdnij_parcial, dnildP, dnivdP, dnildNk, dnivdNk)
        dVtdP, dVtdNk = self.dVt_dP_dNk(kprop, dnldP, dnvdP, dnldNk, dnvdNk, dZldP, dZvdP, dZldNk, dZvdNk, P, Zl, Zv, Nl, Nv)
        import pdb; pdb.set_trace()
        ''' resolução para dVtdNk e dVtdP'''
        #dVt_dNk = self.get_dVt_dNk_analytically(P, Vt, So, l, Nk)
        #No = fprop.phase_mole_numbers[0,0,:]
        #dVt_dP = self.get_dVt_dP_analytically(P, Vt, No, l)
        return dVt_dNk, dVt_dP

    def dZ_dP_dNk(self, dZldP, dZvdP, dZldnij, dZvdnij, dnildP, dnivdP, dnildNk, dnivdNk):
        dZldP = dZldP + np.sum(dZldnij.sum(axis=0) * dnildP, axis = 0)
        dZvdP = dZvdP + np.sum(dZvdnij.sum(axis=0) * dnivdP, axis = 0)
        dZldNk = np.sum(dZldnij[0][:,np.newaxis,:] * dnildNk, axis = 0)
        dZvdNk = np.sum(dZvdnij[0][:,np.newaxis,:] * dnivdNk, axis = 0)
        return dZldP, dZvdP, dZldNk, dZvdNk

    def dVt_dP_dNk(self, kprop, dnldP, dnvdP, dnldNk, dnvdNk, dZldP, dZvdP, dZldNk, dZvdNk, P, Zl, Zv, Nl, Nv):
        coef = (kprop.R * kprop.T / P)
        dVtdP = coef * (Nl * dZldP + Nv * dZvdP + Zl * dnldP + Zv * dnvdP - (Zl * Nl / P + Zv * Nv / P))
        dVtdNk = coef * (Nl * dZldNk + Nv * dZvdNk + Zl * dnldNk + Zv * dnvdNk)
        import pdb; pdb.set_trace()
        #rever isso daqui e talvez as equações novamente
        '''plt.figure(1)
        plt.plot(P / 6894.757, -dVtdP * 6894.757)
        plt.grid()
        plt.show()'''
        return dVtdP, dVtdNk


    def get_phase_derivatives(self, kprop, P, xij, Nj, ph):
        A, B = self.coefficients_cubic_EOS_all(kprop, xij, P)
        Zfunc = np.vectorize(PengRobinson.Z)
        Z = Zfunc(B, A, ph)
        dAdP = self.dA_dP(kprop)
        dBdP = self.dB_dP(kprop)
        dbdnij = self.db_dnij(kprop, Nj, P)
        dadnij = self.da_dnij(kprop, Nj, xij, P)
        dAdnij = self.dA_dnij(kprop, Nj, dadnij, P)
        dBdnij = self.dB_dnij(kprop, Nj, dbdnij, P)
        dZdP_parcial = self.dZ_dP_parcial(dAdP, dBdP, Z, A, B)
        dZdnij_parcial = self.dZ_dnij_parcial(dAdnij, dBdnij, Z, A, B)
        dlnphidP = self.dlnphi_dP(kprop, dAdP, dBdP, dZdP_parcial, Z, A, B)
        dlnphidnij = self.dlnphi_dnij(kprop, dAdnij, dBdnij, dZdnij_parcial, Z, A, B, dbdnij, dadnij, Nj, P)
        dlnfijdP = self.dlnfij_dP(P, dlnphidP)
        dlnfijdnij = self.dlnfij_dnij(kprop, dlnphidnij, xij, Nj)
        return dlnfijdP, dlnfijdnij, dZdP_parcial, dZdnij_parcial, Z #e mais coisa ainda também

    def db_dnij(self, kprop, Nj, P): #OK
        dbdnij = np.empty((kprop.Nc,len(P)))
        dbdnij[:,Nj!=0] = 1 / (Nj[Nj!=0]) * (self.b[:,np.newaxis] - self.bm[Nj!=0])
        dbdnij[:,Nj==0] = 0.
        dbdnij = dbdnij[np.newaxis,:,:]
        return dbdnij

    def da_dnij(self, kprop, Nj, xij, P): #OK
        dadnij = np.empty((kprop.Nc,len(P)))
        dadnij[:,Nj!=0] = 2 / Nj[Nj!=0] * ((xij.T[Nj!=0] @ self.aalpha_ik).T - self.aalpha[Nj!=0])
        dadnij[:,Nj==0] = 0.
        dadnij = dadnij[np.newaxis,:,:]
        return dadnij

    def dA_dnij(self, kprop, Nj, dadnij, P): #OK
        dAdnij = P[np.newaxis,np.newaxis,:] / (kprop.R * kprop.T)**2 * dadnij
        return dAdnij

    def dB_dnij(self, kprop, Nj, dbdnij, P): #OK
        dBdnij = P[np.newaxis,np.newaxis,:] / (kprop.R * kprop.T) * dbdnij
        return dBdnij

    def dA_dP(self, kprop): #OK
        dAdP = self.aalpha / (kprop.R * kprop.T) ** 2
        return dAdP

    def dB_dP(self, kprop): #OK
        dBdP = self.bm / (kprop.R * kprop.T)
        return dBdP

    def dZ_dP_parcial(self, dAdP, dBdP, Z, A, B): #OK
        coef0 = -(Z - B)
        coef1 = -(Z ** 2 - (6 * B + 2) * Z - A + 2 * B + 3 * B **2)
        coef2 = 3 * Z ** 2 - 2 * Z * (1 - B) + (A - 3 * B ** 2 - 2 * B)
        dZdP =  (dAdP * coef0 + dBdP * coef1) / coef2
        return dZdP

    def dZ_dnij_parcial(self,dAdnij, dBdnij, Z, A, B):
        coef0 = -(Z - B)
        coef1 = -(Z ** 2 - (6 * B + 2) * Z - A + 2 * B + 3 * B **2)
        coef2 = 3 * Z ** 2 - 2 * Z * (1 - B) + (A - 3 * B ** 2 - 2 * B)
        dZdnij = (dAdnij * coef0 + dBdnij * coef1) / coef2
        return dZdnij

    def dlnphi_dP(self, kprop, dAdP, dBdP, dZdP, Z, A, B): #OK i guess
        #ai = np.sum(self.aalpha_ik, axis=1)
        coef0 = self.b[:,np.newaxis] / self.bm #coef8
        coef1 = 1 / (Z - B) #coef5
        coef2 = 2 * self.psi / self.aalpha - coef0 #coef0
        coef3 = np.log((Z + (1 + 2**(1/2))*B)/(Z + (1 - 2**(1/2))*B)) #coef1
        coef4 = 2 * (2 **(1/2)) / (Z**2 + 2*Z*B - B**2) #coef2

        dlnphiijdP = coef0 * dZdP - coef1 * (dZdP - dBdP) - coef2 / (2*(2 **(1/2))) * ((B * dAdP - A * dBdP) /  (B **2) * coef3 +
                    A / B *  coef4 * (Z * dBdP - B * dZdP))
        return dlnphiijdP

    def dlnphi_dnij(self, kprop, dAdnij, dBdnij, dZdnij, Z, A, B, dbdnij, dadnij, Nj, P):
        coef0 = self.b[:,np.newaxis] / self.bm #coef8
        coef1 = 1 / (Z - B) #coef5
        coef2 = 2 * self.psi / self.aalpha - coef0 #coef0
        coef3 = np.log((Z + (1 + 2**(1/2))*B)/(Z + (1 - 2**(1/2))*B)) #coef1
        coef4 = 2 * (2 **(1/2)) / (Z**2 + 2*Z*B - B**2) #coef2

        dlnphiijdnij = np.empty((kprop.Nc, kprop.Nc, len(P)))
        dlnphiijdnij[:,:,Nj==0] = 0

        aux1 = coef0[:,np.newaxis,:] / self.bm * (self.bm * dZdnij - (Z - 1) * dbdnij) - coef1 * (dZdnij - dBdnij)
        aux2_1 = coef2[:,np.newaxis,:] * ((B * dAdnij - A * dBdnij) / B**2 * coef3 + A / B * coef4 * (Z * dBdnij - B * dZdnij))
        aux2_2_1 =  2 * (self.aalpha_ik[:,:,np.newaxis] - (Nj / self.aalpha * dadnij + 1) * self.psi[:,np.newaxis,:])
        aux2_2_2 = coef0[:,np.newaxis,:] / self.bm * dbdnij
        dlnphiijdnij[:,:,Nj!=0] = aux1[:,:,Nj!=0] - 1 / (2 * (2 **(1/2))) * (aux2_1[:,:,Nj!=0] + A[Nj!=0] / B[Nj!=0] * coef3[Nj!=0] *
                                (aux2_2_1 [:,:,Nj!=0] / (Nj[Nj!=0] * self.aalpha[Nj!=0]) + aux2_2_2[:,:,Nj!=0]))
        return dlnphiijdnij

    def dlnfij_dP(self, P, dlnphidP):
        dlnfijdP = dlnphidP + 1 / P[np.newaxis,:]
        return dlnfijdP

    def dlnfij_dnij(self, kprop, dlnphidnij, xij, Nj):
        dlnfijdnij = np.empty(dlnphidnij.shape)
        dlnfijdnij[:,:,Nj!=0] = dlnphidnij[:,:,Nj!=0] + 1 / xij[:,np.newaxis,Nj!=0] * (np.identity(kprop.Nc)[:,:,np.newaxis] -
                                xij[:,np.newaxis, Nj!=0]) / Nj[Nj!=0]

        dlnfijdnij[:,:,Nj==0] = 0
        return dlnfijdnij

    def dnij_dNk_dP(self, kprop, dlnfildP, dlnfivdP, dlnfildnij, dlnfivdnij, Nl, Nv):
        dnivdNk = np.empty(dlnfildnij.shape)
        dnivdP = np.empty(dlnfildP.shape)
        dnildNk = np.empty(dlnfildnij.shape)
        dnildP = np.empty(dlnfildP.shape)
        #dlnfivdnij[:,:,Nv==0] = dlnfildnij[:,:,Nv==0]
        #dlnfildnij[:,:,Nl==0] = dlnfivdnij[:,:,Nl==0]

        aux = np.ones(len(Nv),dtype=bool)
        aux[Nv == 0] = False
        aux[Nl == 0] = False

        '''Making dlnf's matrix'''
        dlnfildnij = np.transpose(dlnfildnij[:,:,aux], (2,0,1))
        dlnfivdnij = np.transpose(dlnfivdnij[:,:,aux], (2,0,1))
        matrix = dlnfildnij + dlnfivdnij

        '''Independent tem vector for dnsdP '''
        v = dlnfildP[:,aux] - dlnfivdP[:,aux]

        '''Solving system of equations for dnivdNk and dnivdP'''
        dnivdP_aux = dnivdP[:,aux]
        dnivdNk_aux = dnivdNk[:,:,aux]
        dnivdNk_aux = np.linalg.solve(matrix, dlnfildnij).transpose((1,2,0))
        dnivdP_aux =  np.linalg.solve(matrix, v.T).T

        '''Calculate dnildP and dnildNk'''
        dnivdNk[:,:,aux] = dnivdNk_aux
        dnivdNk[:,:,Nv==0] = 0
        dnivdNk[:,:,Nl==0] = np.identity(kprop.Nc)[:,:,np.newaxis]

        dnivdP[:,aux] = dnivdP_aux
        dnivdP[:,Nv==0] = 0
        dnivdP[:,Nl==0] = 0

        dnildP = - dnivdP
        dnildNk = np.identity(kprop.Nc)[:,:,np.newaxis] - dnivdNk

        dnvdP = np.sum(dnivdP, axis = 0)
        dnldP = - dnvdP
        dnldNk = np.sum(dnildNk, axis = 0)
        dnvdNk = np.sum(dnivdNk, axis = 0)
        return dnldP, dnvdP, dnldNk, dnvdNk, dnildP, dnivdP, dnildNk, dnivdNk

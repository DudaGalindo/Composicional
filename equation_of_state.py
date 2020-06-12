import numpy as np
from math import pi

class PengRobinson:
    def __init__(self, kprop):
        self.coefficientsPR(kprop)

    def coefficientsPR(self, kprop):

        #l - any phase molar composition
        PR_kC7 = np.array([0.379642, 1.48503, 0.1644, 0.016667])
        PR_k = np.array([0.37464, 1.54226, 0.26992])

        k = (PR_kC7[0] + PR_kC7[1] * kprop.w - PR_kC7[2] * kprop.w ** 2 + \
            PR_kC7[3] * kprop.w ** 3) * (1*(kprop.w >= 0.49)) + (PR_k[0] + PR_k[1] * kprop.w - \
            PR_k[2] * kprop.w ** 2) * (1*(kprop.w < 0.49))
        alpha = (1 + k * (1 - (kprop.T / kprop.Tc) ** (1 / 2))) ** 2
        aalpha_i = 0.45724 * (kprop.R * kprop.Tc) ** 2 / kprop.Pc * alpha
        self.b = 0.07780 * kprop.R * kprop.Tc / kprop.Pc
        aalpha_i_reshape = np.ones((kprop.Nc,kprop.Nc)) * aalpha_i[:,np.newaxis]
        self.aalpha_ij = np.sqrt(aalpha_i_reshape.T * aalpha_i[:,np.newaxis]) \
                        * (1 - kprop.Bin)

    def coefficients_cubic_EOS_all(self, kprop, l):
        self.bm = np.sum(l * self.b[:,np.newaxis], axis=0)
        l_reshape = np.ones((self.aalpha_ij).shape)[:,:,np.newaxis] * l[:,np.newaxis,:]
        self.aalpha = (l_reshape * l[np.newaxis,:,:] * self.aalpha_ij[:,:,np.newaxis]).sum(axis=0).sum(axis=0)
        B = self.bm * kprop.P / (kprop.R* kprop.T)
        A = self.aalpha * kprop.P / (kprop.R* kprop.T) ** 2
        self.psi = (l_reshape * self.aalpha_ij[:,:,np.newaxis]).sum(axis = 0).T
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

    def coefficients_cubic_EOS(self, kprop, l):
        self.bm = sum(l * self.b)
        l_reshape = np.ones((self.aalpha_ij).shape) * l[:, np.newaxis]
        self.aalpha = (l_reshape.T * l[:,np.newaxis] * self.aalpha_ij).sum()
        B = self.bm * kprop.P / (kprop.R* kprop.T)
        A = self.aalpha * kprop.P / (kprop.R* kprop.T) ** 2
        self.psi = (l_reshape * self.aalpha_ij).sum(axis = 0)
        return A, B

    def Z(B, A, ph):
        # PR cubic EOS: Z**3 - (1-B)*Z**2 + (A-2*B-3*B**2)*Z-(A*B-B**2-B**3)
        coef = np.array([1, -(1 - B), (A - 2*B - 3*B**2), -(A*B - B**2 - B**3)])

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

    def lnphi(self, kprop, l, ph):
        #l - any phase molar composition
        if l.shape == (len(kprop.P), len(kprop.w)): l = l.T

        A, B = self.coefficients_cubic_EOS_all(kprop,l)
        Zfunc = np.vectorize(PengRobinson.Z)
        Z =Zfunc(B, A, ph)
        #Z = PengRobinson.Z_all(B, A, ph)
        lnphi = self.b / self.bm[:,np.newaxis] * (Z[:,np.newaxis] - 1) - np.log(Z[:,np.newaxis] - B[:,np.newaxis]) - A[:,np.newaxis] / (2 * (2 ** (1/2))
                * B[:,np.newaxis]) * (2 * self.psi / self.aalpha[:,np.newaxis] - self.b / self.bm[:,np.newaxis]) * np.log((Z[:,np.newaxis] + (1 +
                2 ** (1/2)) * B[:,np.newaxis]) / (Z[:, np.newaxis] + (1 - 2 ** (1/2)) * B[:,np.newaxis]))

        return lnphi

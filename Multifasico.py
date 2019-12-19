import numpy as np
from .Checa_Estabilidade import StabilityTest



class Multiphase(object):

    def objective_function_solve2(self,K):
        K1 = max(K);KNc = min(K)

        Lmin = min(0,K1/(K1-1))
        Lmax = max(1,KNc/(KNc-1))
        L = (Lmin+Lmax)/2
        f = 1 #just to get into the loop
        while f >10e-9:

            f = sum((1-K)*self.z/(L+(1-L)*K))
            df = -sum((1-K)**2*self.z/(L+(1-L)*K)**2)

            L = L -f/df
            if L>Lmax: L = Lmax
            elif L<Lmin: L = Lmin

        x = self.z/(L+(1-L)*K)
        y = K*x
        return x,y

    def objective_function_solve(self,K):
        K1 = max(K);KNc = min(K)
        z1 = self.z[K==K1]

        z = self.z[self.z!=z1]; z = z[z!=z[K==KNc]]
        Ki = K[K!=K1]; Ki = Ki[Ki!=KNc]

        x1_min = z1*(1-K1)/(K1-KNc)
        x1_max = (1-KNc)/(K1-KNc)
        x1 = (x1_min+x1_max)/2
        while f >10e-9:

            f = 1 + ((K1-KNc)/(KNc-1))*x1 +\
                sum(((Ki-KNc)/(KNc-1))*\
                z*(K1-1)*x1/((Ki-1)*z1+(K1-Ki)*x1))

            df = ((K1-KNc)/(KNc-1)) +\
                sum(((Ki-KNc)/(KNc-1))*z*z1*(K1-1)*\
                (Ki-1)/((Ki-1)*z1+(K1-Ki)*x1)**2)

            x1 = x1 -f/df
            if x1>x1_max: x1 = x1_max
            elif x1<x1_min: x1 = x1_min

        x = np.zeros(self.Nc-1)
        x = (K1-1)*z*x1/((K-1)*z1+(K1-K)*x1)
        x = np.concatenate(x1,x,axis=None)
        x[self.Nc-1] = 1 - sum(x)
        y = K*x
        return x,y

    def molar_fractions(self):
        # Aqui a correlação de wilson é utilizada apenas para achar o K inicial,
        K = self.K
        fv = 2*np.ones(len(self.z));fl = np.ones(len(self.z))
        while max(abs(fv/fl - 1)) > 10e-9:
            x,y = Multiphase.objective_function_solve2(self,K)
            lnphil = StabilityTest.lnphi(self,x,1)
            lnphiv = StabilityTest.lnphi(self,y,0)
            fv = np.exp(lnphiv)/(y*self.P)
            fl = np.exp(lnphil)/(x*self.P)
            K = (fl/fv)*K
            print(K )
        self.K = K

        print('x',x)
        print('y',y)

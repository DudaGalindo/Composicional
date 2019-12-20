import numpy as np
from .Checa_Estabilidade import StabilityTest



class Multiphase(object):

    def objective_function_solve2(self,K):

        Lmin = min(0,min(K)/(min(K)-1))
        Lmax = max(1,max(K)/(max(K)-1))
        L = (Lmin+Lmax)/2

        f = 1 #just to get into the loop
        i=0
        while abs(f) >10e-9:
            i = i+1
            f = sum((1-K)*self.z/(L+(1-L)*K))
            df = -sum((1-K)**2*self.z/(L+(1-L)*K)**2)
            L = L -f/df
            if L>Lmax: L = Lmax
            elif L<Lmin: L = Lmin

        print(i)
        x = self.z/(L+(1-L)*K)
        y = K*x
        return x,y

    def objective_function_solve(self,K):
        K1 = max(K);KNc = min(K)
        z1 = self.z[K==K1]

        z = self.z[self.z!=z1]
        aux = z[z == self.z[K==KNc]]
        index = np.argwhere(z==aux)
        z = np.delete(z,index[len(index)-1])

        Ki = K[(K!=K1) & (K!=KNc)]

        x1_min = z1*(1-KNc)/(K1-KNc)
        x1_max = (1-KNc)/(K1-KNc)
        x1 = (x1_min+x1_max)/2
        i = 0
        f = 1 #just to get into the first loop
        while abs(f) >10e-9:
            i = i+1
            f = 1 + ((K1-KNc)/(KNc-1))*x1 +\
                sum(((Ki-KNc)/(KNc-1))*\
                z*(K1-1)*x1/((Ki-1)*z1+(K1-Ki)*x1))

            df = ((K1-KNc)/(KNc-1)) +\
                sum(((Ki-KNc)/(KNc-1))*z*z1*(K1-1)*\
                (Ki-1)/((Ki-1)*z1+(K1-Ki)*x1)**2)

            x1 = x1 -f/df
            if (x1>x1_max) | (x1<x1_min): x1 = (x1_min+x1_max)/2
        print(i)

        xi = (K1-1)*z*x1/((Ki-1)*z1+(K1-Ki)*x1)
        xNc = 1 - sum(xi) - x1
        x = np.concatenate((x1,xi,xNc),axis=None)
        y = K*x
        return x,y

    def molar_fractions(self):
        # Aqui a correlação de wilson é utilizada apenas para achar o K inicial,
        K = self.K
        fv = 2*np.ones(len(self.z));fl = np.ones(len(self.z))
        while max(abs(fv/fl - 1)) > 10e-9:
            x,y = Multiphase.objective_function_solve(self,K)
            lnphil = StabilityTest.lnphi(self,x,1)
            lnphiv = StabilityTest.lnphi(self,y,0)
            fv = np.exp(lnphiv)*(y*self.P)
            fl = np.exp(lnphil)*(x*self.P)
            K = (fl/fv)*K
        print(x,y)
        self.K = K

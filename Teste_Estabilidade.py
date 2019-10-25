import numpy as np
import math
from math import pi
import thermo
import matplotlib.pyplot as plt

## Encontrar os pontos estacionarios. Estes correspondem aos pontos nos quais a derivada de g com respeito a Y é 0
class fugacity:
    def coefficientsPR(w,Bin,R,Tc,Pc,T,P,Nc):
        a = np.zeros(Nc)
        b = np.zeros(Nc)
        K = np.zeros(Nc)
        aalpha_i = np.zeros(Nc)
        aalpha_ij = np.zeros((Nc,Nc))

        for i in range(0,Nc):
            k = 0.37464 + 1.5422*w[i] -0.26992*w[i]**2; #kappa - PR eos
            alpha = (1+k*(1-(T/Tc[i])**(1/2)))**2; #coef - PR eos
            aalpha_i[i] = 0.45724*R**2*Tc[i]**2*alpha/Pc[i]; #coef - PR eos
            b[i] = 0.07780*R*Tc[i]/Pc[i]; #coef - PR eos
            K[i] = math.exp(5.37*(1+w[i])*(1-1/(T/Tc[i])))/(P/Pc[i]); #acho que isso encontro em algum livro
        return K,b,aalpha_ij,aalpha_i



    def V_PR(bm,R,T,P,A,aalpha,ph):
        ## Acredito que aqui ele tenha usado a EOS de PR da equação cúbica pra resolver, identificando se vai ter uma,duas ou tres raízes.
        # Não sei de onde ele tirou essa forma de escrever, mas pretendo testar com afunção do thermo que já retorna Vl e Vv

        p = bm - R*T/P
        q = -3*bm**2 - 2*R*T*bm/P + aalpha/P
        r = bm**3 + R*T*bm**2/P - aalpha*bm/P
        a1 = 1/3*(3*q-p**2)
        a2 = 1/27*(2*p**3 - 9*p*q + 27*r)

        TEST = (a2**2/4 + a1**3/27)
        if TEST < 0:
           phi = math.acos(-(a2/2)/(-a1**3/27)**0.5)
           V = np.zeros(3)
           V[0] = 2*(-a1/3)**0.5*math.cos(phi/3) - p/3
           V[1] = 2*(-a1/3)**0.5*math.cos(phi/3 + 2*pi/3) -p/3
           V[2] = 2*(-a1/3)**0.5*math.cos(phi/3 + 4*pi/3) -p/3
           Vl = min(V)
           Vv = max(V)
           roots = 3
        elif TEST == 0:
            V = np.zeros(2)
            b1 = (-a2/2 + (a2**2/4 + a1**3/27)**(1/2))**(1/3)
            V[0] = 2*b1 - p/3
            V[1] = -b1 - p/3
            Vl = min(V)
            Vv = max(V)
            roots=2
        else:
            b1 = (-a2/2 + (a2**2/4 + a1**3/27)**(1/2))**(1/3)
            b22 =  -a2/2 - (a2**2/4 + a1**3/27)**(1/2)
            if b22 < 0: aux = -(-b22)*(1/3)
            else: aux = b22**(1/3)
            b2 = np.sign(( -a2/2 - (a2**2/4 + a1**3/27)**(1/2)))*abs(aux)
            if ph==1:
               Vl = b1 + b2 - p/3
               Vv = 1e-20
            else:
                Vv = b1 + b2 - p/3
                Vl = 1e-20
            roots=1

        return Vl,Vv


    def lnphi(x,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_i,aalpha_ij,b,ph): #ph stands for phase (1-liquid, 2-vapor)
        lnphi = np.zeros(Nc)
        lnphi_termo = np.zeros(Nc)

        aalpha = 0; bm = 0
        for i in range(0,Nc):
            bm = bm + x[i]*b[i]
            for j in range(0,Nc):
                aalpha_ij[i,j] = (aalpha_i[i]*aalpha_i[j])**(1/2) * (1-Bin[i,j]); # aalphaij
                aalpha = aalpha + x[i]*x[j]*aalpha_ij[i,j]
        B = bm*P/(R*T)
        A = aalpha*P/(R*T)**2

        #Teste do cálculo do volume da fase:
        '''CHECAR SE O MÉTODO DO THERMO TA BATENDO COM O QUE O DR. Schmall ESCREVEU'''
        eos = thermo.eos_mix.PRMIX(Tc,Pc,w,x,Bin,T=T,P=P)
        Vl_t = eos.V_l
        #Vv_t = eos.V_g

        Vl,Vv = fugacity.V_PR(bm,R,T,P,A,aalpha,ph)

        Zv = P*Vv/(R*T)
        Zl = P*Vl/(R*T)
        for i in range(0,Nc):
            psi_i = 0
            for j in range(0,Nc):
                psi_i = psi_i + x[j]*aalpha_ij[i,j]
            if ph == 1:
                lnphi[i] = b[i]/bm*(Zl-1) - np.log(Zl - B) - \
                    A/(2*(2**(1/2))*B)*(2*psi_i/aalpha -\
                    b[i]/bm)*np.log((Zl+(1+2**(1/2))*B)/(Zl+(1-2**(1/2))*B))
            else:
                    lnphi[i] = b[i]/bm*(Zv-1) - np.log(Zv - bm*P/(R*T)) -\
                    A/(2*(2**(1/2))*B)*(2*psi_i/aalpha -\
                    b[i]/bm)*np.log((Zv+(1+2**(1/2))*B)/(Zv+(1-2**(1/2))*B))
        return lnphi



class StabilityTest:
    def Stability(w,Bin,R,Tc,Pc,T,P,Nc,z):
        K,b,aalpha_ij,aalpha_i = fugacity.coefficientsPR(w,Bin,R,Tc,Pc,T,P,Nc)
        #****************************Initial guess******************************#
        ## Both approaches bellow should be used in case the phase is in the critical region
        ''' VER QUAL CRITÉRIO RETORNA SE É ESTÁVEL OU NÃO'''
        ## Used alone when the phase investigated is liquid like (ph == 1)
        Y = z/K
        Y_old = 0.9*Y #just to pass the first while below
        y = Y/sum(Y)

        while max(abs(Y/Y_old - 1))>1e-9: #convergência - Y_old/Y termo a termo
            Y_old = Y
            lnphiz = fugacity.lnphi(z,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_i,aalpha_ij,b,1); #by the PR EOS
            lnphiy = fugacity.lnphi(y,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_i,aalpha_ij,b,2);
            for i in range(len(Y)):
                Y[i] = math.exp( math.log(z[i]) + lnphiz[i] - lnphiy[i] )
            y = Y/sum(Y);

        if sum(Y) <= 1: return 'estavel'
        else: return'instavel'

        ## Used alone when the phase investigated is vapour like (ph == 2)
        Y = K*z
        Y_old = 0.9*Y
        y = np.divide(Y,sum(Y)) #ta dando erro

        while max(abs(Y/Y_old - 1))>1e-9:
            Y_old = Y
            lnphiz = fugacity.lnphi(z,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_i,aalpha_ij,b,2);
            lnphiy = fugacity.lnphi(y,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_i,aalpha_ij,b,1);
            for i in range(len(Y)):
                Y[i] = math.exp( math.log(z[i]) + lnphiz[i] - lnphiy[i] )
            y = np.divide(Y,sum(Y))

        if sum(Y) <= 1: return 'estavel'
        else: return'instavel'


    def TPD(Nc,T,P,R,Tc,Pc,Bin,w):
        #*************************Tangent Plane distance plot******************#
        t = np.linspace(0.01,0.99,0.9/0.02) #vetor auxiliar
        print(t)
        TPD = np.zeros(len(t)) ##F
        for i in range(0,len(t)):
            aux = 0;
            lnphiz = StabilityTest.lnphi(z,Nc,T,P,R,Tc,Pc,Bin,w,1) #original phase
            x = np.array([1-t[i],t[i]]) #new phase composition (1-t e t) - apenas válido para Nc=2 acredito eu
            lnphix = StabilityTest.lnphi(x,Nc,T,P,R,Tc,Pc,Bin,w,2) #new phase (vapor- ph=2)
            for j in range(0,Nc):
                fix = math.exp(lnphix[j])*x[j]*P
                fiz = math.exp(lnphiz[j])*z[j]*P; #se acaba dividindo os dois embaixo, não tem pra quê multiplicar pela pressão -ps: realmente nem sei o que a pressão ta fazendo aí
                aux = aux + x[j]*R*T*(math.log(fix/fiz))
            TPD[i] = aux

        plt.figure(0)
        plt.plot(t,TPD)
        plt.xlabel('x')
        plt.ylabel('TPD')
        plt.show()

        return TPD

import numpy as np
import math
from math import pi
import thermo
import matplotlib.pyplot as plt

## Encontrar os pontos estacionarios. Estes correspondem aos pontos nos quais a derivada de g com respeito a Y é 0
## Todas as equações foram confirmadas pelo livro do Dandekar e a biblioteca do thermo

class fugacity:
    def coefficientsPR(w,Bin,R,Tc,Pc,T,P,Nc):
        a = np.zeros(Nc)
        b = np.zeros(Nc)
        K = np.zeros(Nc)
        aalpha_i = np.zeros(Nc); aalpha_ij= np.zeros((Nc,Nc))


        for i in range(0,Nc):
            k = 0.37464 + 1.54226*w[i] -0.26992*w[i]**2;
            alpha = (1+k*(1-(T/Tc[i])**(1/2)))**2;
            aalpha_i[i] = 0.45724*(R*Tc[i])**2*alpha/Pc[i];
            b[i] = 0.07780*R*Tc[i]/Pc[i];
            K[i] = math.exp(5.37*(1+w[i])*(1-1/(T/Tc[i])))/(P/Pc[i]); # Wilson equation (K - equilimbrium ratio)

        for i in range(0,Nc):
            for j in range(0,Nc):
                aalpha_ij[i,j] = (aalpha_i[i]*aalpha_i[j])**(1/2) *(1.
                                  - Bin[i,j]); # aalphaij
        return K,b,aalpha_ij



    '''def V_PR(bm,R,T,P,A,aalpha,ph):
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
            b2 = np.sign(( -a2/2 - (a2**2/4 + a1**3/27)**(1/2)))*abs(-a2/2 - (a2**2/4 + a1**3/27)**(1/2))**(1/3)

            if ph==1:
               Vl = b1 + b2 - p/3
               Vv = 1e-20
            else:
                Vv = b1 + b2 - p/3
                Vl = 1e-20
            roots=1

        return Vl,Vv'''


    def lnphi(x,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b):
        lnphi = np.zeros(Nc)
        lnphi_termo = np.zeros(Nc)
        aalpha = 0; bm = 0;

        bm = sum(x*b)
        B = bm*P/(R*T)
        for i in range(0,Nc):
            for j in range(0,Nc):
                aalpha = aalpha + x[i]*x[j]*aalpha_ij[i,j]
        A = aalpha*P/(R*T)**2

        eos = thermo.eos_mix.PRMIX(Tcs=Tc,Pcs=Pc,omegas=w,zs=x,kijs=Bin,T=T,P=P)
        if eos.phase == 'l': V = eos.V_l
        else: V = eos.V_g
        
        Z = P*V/(R*T)

        for i in range(0,Nc):
            psi_i = 0
            for j in range(0,Nc):
                psi_i = psi_i + x[j]*aalpha_ij[i,j]

            lnphi[i] = b[i]/bm*(Z-1) - np.log(Z - B) - \
                A/(2*2**(1/2)*B)*(2*psi_i/aalpha -\
                b[i]/bm)*np.log((Z+(1+2**(1/2))*B)/(Z+(1-2**(1/2))*B))

        return lnphi



class StabilityTest:
    def Stability(w,Bin,R,Tc,Pc,T,P,Nc,z):
        K,b,aalpha_ij = fugacity.coefficientsPR(w,Bin,R,Tc,Pc,T,P,Nc)

    #****************************Initial guess******************************#
    ## Both approaches bellow should be used in case the phase is in the critical region

        #identifying the initial phase 'z' mole fraction
        eos = thermo.eos_mix.PRMIX(Tcs=Tc,Pcs=Pc,omegas=w,zs=z,kijs=Bin,T=T,P=P)

        lnphiz = fugacity.lnphi(z,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b)
        ## Used alone when the phase investigated is clearly liquid like (ph == 1)
        if eos.phase == 'l' or eos.phase == 'l/g':
            Y = z/K
            Y_old = 0.9*Y
            y = Y/sum(Y)

            while max(abs(Y/Y_old - 1))>1e-10: #convergência - critério do Schmall - vi no livro do Dandekar que tinha outro que segundo ele era melhor
                Y_old = Y
                lnphiy = fugacity.lnphi(y,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b)

                for i in range(len(Y)):
                    Y[i] = math.exp( math.log(z[i]) + lnphiz[i] - lnphiy[i] )

                y = Y/sum(Y);

            if sum(Y) <= 1: print('estavel')#estavel2 = 1
            else: print('instavel') #estavel2 = 0

    ## Used alone when the phase investigated is clearly vapour like (ph == 2)
        if eos.phase == 'g' or eos.phase == 'l/g':
            Y = K*z
            Y_old = 0.9*Y
            y = np.divide(Y,sum(Y)) #ta dando erro

            while max(abs(Y/Y_old - 1))>1e-9:
                Y_old = Y
                lnphiy = fugacity.lnphi(y,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b);
                for i in range(len(Y)):
                    Y[i] = math.exp( math.log(z[i]) + lnphiz[i] - lnphiy[i] )
                y = np.divide(Y,sum(Y))

            if sum(Y) <= 1: print('estavel')#estavel2 = 1
            else: print('instavel') #estavel2 = 0

        '''
        if estavel1==1 and estavel2==1: return 'estavel'
        else: return 'instavel'''

    def TPD(Nc,T,P,R,Tc,Pc,Bin,w,z): #atualmente só funciona para 2D
        x = np.zeros()
        K,b,aalpha_ij = fugacity.coefficientsPR(w,Bin,R,Tc,Pc,T,P,Nc)
        #**********************Tangent Plane distance plot*********************#
        t = np.linspace(0.01,0.99,0.9/0.002) #vetor auxiliar
        TPD = np.zeros(len(t)) ##F
        for i in range(0,len(t)):
            aux = 0;
            lnphiz = fugacity.lnphi(z,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b) #original phase

            x = np.array([1-t[i],t[i]]) #new phase composition (1-t e t) - apenas válido para Nc=2 acredito eu.

            '''O modo que x varia implica no formato de TPD. No presente exemplo,
            a fração molar do segundo componente de x varia direto com t, que é a
            variável de plotagem. Logo, a distancia dos planos tangentes será
            zeros em z[Nc-1]. O contrário ocorreria'''
            lnphix = fugacity.lnphi(x,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b); #new phase (vapor- ph=2)
            for j in range(0,Nc):
                fix = math.exp(lnphix[j])*x[j]*P
                fiz = math.exp(lnphiz[j])*z[j]*P
                aux = aux + x[j]*R*T*(math.log(fix/fiz))
            TPD[i] = aux

        for i in range(len(t)):
            if abs(TPD[i] - K[1]) < 0.1:
                print(TPD[i])
                print('ysp',t[i])

        plt.figure(0)
        plt.plot(t,TPD)
        plt.xlabel('x')
        plt.ylabel('TPD')
        plt.show()

        return TPD

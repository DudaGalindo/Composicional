import numpy as np
import math
from math import pi
import thermo
import matplotlib.pyplot as plt

## Encontrar os pontos estacionarios. Estes correspondem aos pontos nos quais a derivada de g com respeito a Y é 0
## Todas as equações foram confirmadas pelo livro do Dandekar e a biblioteca do thermo

class fugacity:
    def coefficientsPR(w,Bin,R,Tc,Pc,T,P,Nc,C7):
        a = np.zeros(Nc)
        b = np.zeros(Nc)
        K = np.zeros(Nc)
        aalpha_i = np.zeros(Nc); aalpha_ij= np.zeros((Nc,Nc))


        for i in range(0,Nc):
            if C7 == 'y':# For heavier components:
                k = 0.379642+1.48503*w[i]-0.1644*w[i]**2+0.016667*w[i]**3
            else: k = 0.3746 + 1.54226*w[i] -0.26992*w[i]**2;
            alpha = (1+k*(1-(T/Tc[i])**(1/2)))**2;
            aalpha_i[i] = 0.45724*(R*Tc[i])**2*alpha/Pc[i];
            b[i] = 0.07780*R*Tc[i]/Pc[i];
            K[i] = math.exp(5.37*(1+w[i])*(1-1/(T/Tc[i])))/(P/Pc[i]); # Wilson equation (K - equilimbrium ratio)

        for i in range(0,Nc):
            for j in range(0,Nc):
                aalpha_ij[i,j] = (aalpha_i[i]*aalpha_i[j])**(1/2) *(1.
                                  - Bin[i,j]); # aalphaij
        return K,b,aalpha_ij


    def Z_PR(B,A,ph):
        #Z**3 - (1-B)*Z**2 + (A-2*B-3*B**2)*Z-(A*B-B**2-B**3)
        coef = [1,-(1-B),(A-2*B-3*B**2),-(A*B-B**2-B**3)]
        Z = np.roots(coef)
        root = np.isreal(Z) # return True for real roots

        # Saving them in the result vector Zr
        pos = np.where(root == True) #position where the real roots are
        Z_reais = np.zeros(len(pos))

        for i in range(0,len(Z_reais)):
            Z_reais[i] = np.real(Z[pos[i]]) #Saving the real roots

        if sum(pos)>1:
            if ph == 'l': Z_ans = min(Z_reais)
            else: Z_ans = max(Z_reais)
        else: Z_ans = Z_reais

        '''OBS: O caso acima considera que pode haver mais de uma raiz real uma
        vez que o dado de entrada pode ser uma substancia simples.'''

        return Z_ans

    def V_PR(bm,R,T,P,A,aalpha,ph):
        ## Acredito que aqui ele tenha usado a EOS de PR da equação cúbica pra resolver, identificando se vai ter uma,duas ou tres raízes.
        # Não sei de onde ele tirou essa forma de escrever

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

            if ph=='l':
               Vl = b1 + b2 - p/3
               Vv = 1e-20
            else:
                Vv = b1 + b2 - p/3
                Vl = 1e-20
            roots=1

        return Vl,Vv


    def lnphi(x,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b,ph):
        lnphi = np.zeros(Nc)
        lnphi_termo = np.zeros(Nc)
        aalpha = 0; bm = 0;

        bm = sum(x*b)
        B = bm*P/(R*T)
        for i in range(0,Nc):
            for j in range(0,Nc):
                aalpha = aalpha + x[i]*x[j]*aalpha_ij[i,j]
        A = aalpha*P/(R*T)**2

        '''eos = thermo.eos_mix.PRMIX(Tcs=Tc,Pcs=Pc,omegas=w,zs=x,kijs=Bin,T=T,P=P)
        if eos.phase == 'l': V = eos.V_l
        else: V = eos.V_g'''
        '''Vl,Vv = fugacity.V_PR(bm,R,T,P,A,aalpha,ph)
        if ph == 'l': V = Vl
        else: V = Vv
        Z = P*V/(R*T)'''
        Z = fugacity.Z_PR(B,A,ph)

        for i in range(0,Nc):
            psi_i = 0
            for j in range(0,Nc):
                psi_i = psi_i + x[j]*aalpha_ij[i,j]

            lnphi[i] = b[i]/bm*(Z-1) - np.log(Z - B) - \
                A/(2*2**(1/2)*B)*(2*psi_i/aalpha -\
                b[i]/bm)*np.log((Z+(1+2**(1/2))*B)/(Z+(1-2**(1/2))*B))

        return lnphi



class StabilityTest:
    def Stability(w,Bin,R,Tc,Pc,T,P,Nc,C7,z):
        K,b,aalpha_ij = fugacity.coefficientsPR(w,Bin,R,Tc,Pc,T,P,Nc,C7)

    #****************************Initial guess******************************#
    ## Both approaches bellow should be used in case the phase is in the critical region

        #identifying the initial phase 'z' mole fraction
        '''eos = thermo.eos_mix.PRMIX(Tcs=Tc,Pcs=Pc,omegas=w,zs=z,kijs=Bin,T=T,P=P)
        print(eos.phase)'''
        ''' Talvez aqui eu possa fazer uma modificação no que diz respeito a: o
        estado inicial, se eu conheço ele, sei qual sua fase (liq ou vapor),
        ja me limita na escolha dos testes - pensar mais sobre isso'''

        ## Used alone when the phase investigated is clearly vapor like (ph == g)

        Y = z/K
        Yold = 0.9*Y
        y = Y/sum(Y)

        while max(abs(Y/Yold - 1))>1e-10: #convergência
            Yold = np.copy(Y)
            lnphiz = fugacity.lnphi(z,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b,'l')
            lnphiy = fugacity.lnphi(y,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b,'g')

            for i in range(Nc):
                Y[i] = math.exp( math.log(z[i]) + lnphiz[i] - lnphiy[i] )
            y = Y/sum(Y);
        print(sum(Y))
        if sum(Y) <= 1: print('estavel')
        else: print('instavel') #estavel2 = 0

        ''' If this first test returns stable: The original mixture corresponds
        that the Gibbs free energy is at global minimum, there is only one
        stationary point is where y = z. So, it does not have a new phase,
        doesn't need a phase split. If it returns unstable: There is another
        composition that makes the Gibbs free energy at his global minimum wich
        corresponds to a negative value where sum(Y)>1'''

        ## Used alone when the phase investigated is clearly liquid like (ph == l)

        Y = K*z
        Y_old = 0.9*Y
        y = Y/sum(Y)

        while max(abs(Y/Y_old - 1))>1e-10:
            Y_old = np.copy(Y)
            lnphiz = fugacity.lnphi(z,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b,'g')
            lnphiy = fugacity.lnphi(y,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b,'l')
            
            for i in range(len(Y)):
                Y[i] = math.exp( math.log(z[i]) + lnphiz[i] - lnphiy[i] )
            y = Y/sum(Y)
        print(sum(Y))
        if sum(Y) <= 1: print('estavel')#estavel2 = 1
        else: print('instavel') #estavel2 = 0

    ''''The same thing happens here. The difference is that, the original phase
    is gas, and then the "new" phase is liquid. ''' #completar aqui ainda.

    def TPD(Nc,C7,T,P,R,Tc,Pc,Bin,w,z): #atualmente só funciona para 2D
        x = np.zeros(Nc)
        K,b,aalpha_ij = fugacity.coefficientsPR(w,Bin,R,Tc,Pc,T,P,Nc,C7)
        #**********************Tangent Plane distance plot*********************#
        t = np.linspace(0.01,0.99,0.9/0.002) #vetor auxiliar
        TPD = np.zeros(len(t)) ##F
        for i in range(0,len(t)):
            aux = 0;
            lnphiz = fugacity.lnphi(z,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b,'l') #original phase

            x = np.array([1-t[i],t[i]]) #new phase composition (1-t e t) - apenas válido para Nc=2 acredito eu.

            '''O modo que x varia implica no formato de TPD. No presente exemplo,
            a fração molar do segundo componente de x varia direto com t, que é a
            variável de plotagem. Logo, a distancia dos planos tangentes será
            zeros em z[Nc-1]. O contrário ocorreria'''
            lnphix = fugacity.lnphi(x,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b,'g'); #new phase (vapor- ph=2)
            for j in range(0,Nc):
                fix = math.exp(lnphix[j])*x[j]*P
                fiz = math.exp(lnphiz[j])*z[j]*P
                aux = aux + x[j]*R*T*(math.log(fix/fiz))
            TPD[i] = aux

        '''for i in range(len(t)):
            if (TPD[i]) < 0:
                print(TPD[i])
                print('ysp',t[i])'''

        plt.figure(0)
        plt.plot(t,TPD)
        plt.xlabel('x')
        plt.ylabel('TPD')
        plt.show()

        return TPD

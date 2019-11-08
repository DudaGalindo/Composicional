import numpy as np
import math
from math import pi
import thermo
import matplotlib.pyplot as plt

## Encontrar os pontos estacionarios. Estes correspondem aos pontos nos quais a derivada de g com respeito a Y é 0
## Todas as equações foram confirmadas pelo livro do Dandekar e a biblioteca do thermo
'''Some observations:
    This code only works for phase by phase (one phase at a time). For example,
    if the mixture is in the two phase region, it should be identified the molar
    fractions of the components in each phase.'''

class Fugacity:
    def coefficientsPR(w,Bin,R,Tc,Pc,T,P,Nc,C7):
        # I think that the method to calculate k would be just the general one for heavier components
        k = (0.379642+1.48503*w-0.1644*w**2+0.016667*w**3)*C7+\
            (0.3746 + 1.54226*w -0.26992*w**2)*(1-C7)
        alpha = (1+k*(1-np.sqrt(T/Tc)))**2;
        aalpha_i = 0.45724*(R*Tc)**2/Pc*alpha
        b = 0.07780*R*Tc/Pc;
        K = np.exp(5.37*(1+w)*(1-Tc/T))*(Pc/P); # Wilson equation (K - equilimbrium ratio)
        aalpha_i_reshape = np.ones((Nc,Nc))*aalpha_i[:,np.newaxis]
        aalpha_ij = np.sqrt(aalpha_i_reshape.T*aalpha_i[:,np.newaxis])*(1.- Bin)

        return K,b,aalpha_ij


    def Z_PR(B,A,ph):
        # PR cubic EOS: Z**3 - (1-B)*Z**2 + (A-2*B-3*B**2)*Z-(A*B-B**2-B**3)
        coef = [1,-(1-B),(A-2*B-3*B**2),-(A*B-B**2-B**3)]
        Z = np.roots(coef)
        root = np.isreal(Z) # return True for real roots
        real_roots_position = np.where(root == True) #position where the real roots are - criada apenas para organização
        Z_reais = np.real(Z[real_roots_position[:]]) #Saving the real roots

        ''' This part below, considers that the phase is composed by a pure
         component, so the EOS model can return more than one real root.
            If liquid, Zl = min(Z) and gas Zv = max(Z).
            You can notice that, if there's only one real root, it works as well.'''

        Z_ans = min(Z_reais)*ph + max(Z_reais)*(1-ph)

        return Z_ans


    def lnphi(x,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b,ph):
        bm = sum(x*b)
        B = bm*P/(R*T)
        x_reshape = np.ones(aalpha_ij.shape)*x[:,np.newaxis]
        aalpha = (x_reshape.T*x[:,np.newaxis]*aalpha_ij).sum()
        A = aalpha*P/(R*T)**2

        Z = Fugacity.Z_PR(B,A,ph)

        psi = (x_reshape*aalpha_ij).sum(axis = 0)
        lnphi = b/bm*(Z-1)-np.log(Z-B)-A/(2*(2**(1/2))*B)*\
                (2*psi[:]/aalpha-b/bm)*np.log((Z+(1+2**(1/2))*B)/\
                (Z+(1-2**(1/2))*B))

        return lnphi



class StabilityTest:
    def Stability(w,Bin,R,Tc,Pc,T,P,Nc,C7,z):
        K,b,aalpha_ij = Fugacity.coefficientsPR(w,Bin,R,Tc,Pc,T,P,Nc,C7)

    #****************************INITIAL GUESS******************************#
    ## Both approaches bellow should be used in case the phase is in the critical region

        #identifying the initial phase 'z' mole fraction
        ''' In the lnphi function: 0 stands for gas phase and 1 for liquid'''

    #*****************************Test one**********************************#
        #Used alone when the phase investigated (y) is clearly vapor like (ph == g)

        Y = z/K
        Yold = 0.9*Y
        y = Y/sum(Y)

        while max(abs(Y/Yold - 1))>1e-10: #convergência
            Yold = np.copy(Y)
            lnphiz = Fugacity.lnphi(z,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b,1)
            lnphiy = Fugacity.lnphi(y,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b,0)
            Y = np.exp(np.log(z) + lnphiz - lnphiy)
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
    #*****************************Test two**********************************#
        #Used alone when the phase investigated (y) is clearly liquid like (ph == l)

        Y = K*z
        Y_old = 0.9*Y
        y = Y/sum(Y)

        while max(abs(Y/Y_old - 1))>1e-10:
            Y_old = np.copy(Y)
            lnphiz = Fugacity.lnphi(z,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b,0)
            lnphiy = Fugacity.lnphi(y,Nc,T,P,R,Tc,Pc,Bin,w,aalpha_ij,b,1)
            Y = np.exp(np.log(z) + lnphiz - lnphiy)
            y = Y/sum(Y)
        print(sum(Y))
        if sum(Y) <= 1: print('estavel')#estavel2 = 1
        else: print('instavel') #estavel2 = 0

        '''The same thing happens here. The difference is that, the original
        phase is gas, and then the "new" phase is liquid. '''

    '''OBS: In both tests, the fugacity coefficient of the original phase is
    calculated. For a mixture, however, they have the same value once there
    is only one compressibility factor. The ph = 'g' or ph = 'l' just inffer if
    the phase of the original mixture is unknown, so the initial guess changes,
    or if the present fluid is composed by only one component (its not a mixture).
    This means that, if it's a mixture and you know the phase is liquid for example,
    and in a way ir returns a instability in the second test, you will probably have
    a two phase liquid system. '''#-Discutir com alguém sobre essa última afirmação


    def TPD(Nc,C7,T,P,R,Tc,Pc,Bin,w,z): #atualmente só funciona para 2D - falta modificar ainda. Mas só modifico quando eu ver onde vou usar isso
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

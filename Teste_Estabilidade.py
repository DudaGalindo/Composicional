import numpy as np
import math
from math import pi
import thermo

class fugacity:
    def lnphi(x,a,b,A,Nc,T,P,R,ph): #ph stands for phase (1-liquid, 2-vapor)
        lnphi = np.zeros(Nc)
        lnphi_termo = np.zeros(Nc)
        Am = 0
        Bm = 0
        for i in range(0,Nc):
            for j in range(0,Nc):
                Am = Am + x[i]*x[j]*A[i,j]
            Bm = Bm + x[i]*b[i]
            p = Bm - R*T/P
            q = -3*Bm**2 - 2*R*T*Bm/P + Am/P
            r = Bm**3 + R*T*Bm**2/P - Am*Bm/P
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

            Zv = P*Vv/(R*T)
            Zl = P*Vl/(R*T)

            for i in range(0,Nc):
                if ph ==1:
                    Z = Zl
                else: Z = Zv
                lnphi_thermo[i] = thermo.eos_mix.PRMIX.fugacity_coefficients(Z,zs=x) #PARA COMPARAR COM A FUNÇÃO DO THERMO

                xA = 0
                for j in range(0,Nc):
                    xA = xA + x[j]*A[i,j]

                if ph == 1:
                        lnphi[i] = b[i]/Bm*(Zl-1) - np.log(Zl - Bm*P/(R*T)) - \
                        Am/(2*2**(1/2)*Bm*R*T)*(2*xA/Am -\
                        b[i]/Bm)*np.log((Zl+(1+2**(1/2))*Bm*P/(R*T))/(Zl+(1-2**(1/2))*Bm*P/(R*T)))
                else:
                        lnphi[i] = b[i]/Bm*(Zv-1) - np.log(Zv - Bm*P/(R*T)) -\
                        Am/(2*2**(1/2)*Bm*R*T)*(2*xA/Am -\
                        b[i]/Bm)*np.log((Zv+(1+2**(1/2))*Bm*P/(R*T))/(Zv+(1-2**(1/2))*Bm*P/(R*T)))
        return lnphi,lnphi_thermo





# Dados de entrada:
z = np.array([0.2,0.8]) #exemplo aleatório
caso = 2

if caso == 2:
    Nc = 2;
    P = (100 * 10E5 / 101325) * 14.7;
    R = 10.73159;
    T = 350 *(9/5);
    Tc = (9/5)*np.array([190.6, 460.4]);
    Pc = 14.7*np.array([45.4,33.4]);
    w = np.array([0.008,0.227]);
    Bin = np.array([[0,0.0236],[0.0236,0]]);

a = np.zeros(Nc)
b = np.zeros(Nc)
K = np.zeros(Nc)
A = np.zeros([Nc,Nc])

for i in range(0,Nc):
    k = 0.37464 + 1.5422*w[i] -0.26992*w[i]**2; #kappa - PR eos
    alpha = (1+k*(1-(T/Tc[i])**(1/2)))**2; #coef - PR eos
    a[i] = 0.45724*R**2*Tc[i]**2*alpha/Pc[i]; #coef - PR eos
    b[i] = 0.07780*R*Tc[i]/Pc[i]; #coef - PR eos
    K[i] = math.exp(5.37*(1+w[i])*(1-1/(T/Tc[i])))/(P/Pc[i]); #acho que isso encontro em algum livro

for i in range(0,Nc):
    for j in range(0,Nc):
        A[i,j] = (a[i]*a[j])**(1/2) * (1-Bin[i,j]); #acho que isso encontro em algum livro


## Encontrar os pontos estacionarios. Estes correspondem aos pontos nos quais a derivada de g com respeito a Y é 0

#****************************Estimativas iniciais******************************#
## Both approaches bellow should be used in case the phase is in the critical region
## Used alone when the phase investigated is liquid like
Y = z/K
Y_old = 0.9*Y #só para passar do primeiro while
y = Y/sum(Y)
#lembrete: phi - coeficientes de fugacidade

#while max(abs(Y_old/Y - 1))>1e-9: #convergência - Y_old/Y termo a termo
#    Y_old = Y
lnphiz,lnphiz_thermo = fugacity.lnphi(z,a,b,A,Nc,T,P,R,1); #ver quais os métodos utilizados para tanto
lnphix,lnphix_thermo = fugacity.lnphi(y,a,b,A,Nc,T,P,R,2);
print(lnphix,lnphix_thermo)
print(lnphiz,lnphiz_thermo)
for i in range(len(Y)):
    Y[i] = math.exp(-lnphix[i] + math.log(z[i]) + lnphiz[i])
y = Y/sum(Y)
print('lnphiz,lnphix',lnphiz,lnphix)

## Used alone when the phase investigated is vapour like
#Y = np.matmul(K,z)
#Y_old = 0.9*Y

## Routine for implementing IMPEC compositional formulation

class IMPEC:
    def pressure(n,condicoesDeContorno,permeab_matriz,viscosidade,num_phases):
        return p

    def Darcy(n_phases,p,pc,D,phase_density):
        return phase_velocities

    def composition(propC[k],phase_velocities,num_phases):
        return Nc

    def saturation(num_phases, p, p_old):
        return S

    def times(t,tf,CFL,h,v): ##DÚVIDA
        return deltaT

    def IMPEC(t,tf,n):
         # n = número de elementos da minha malha - malha 1D uniforme e estruturada
         # dados: Nc[k] inicial (de k=0 a k componentes)
         #        propriedade dos componentes (ex: fração molar do componente na fase (xj ou yj onde j = fase),
         #                                         densidade molar da fase
         #                                         saturação inicial das fases
         #                                         initial pressure = p_old

        deltaT = IMPEC.times(t,tf,CFL,h,v)
        t = t + deltaT

        while t < tfinal: #contador de tempo de t até tfinal

            p = IMPEC.pressure(n,condicoesDeContorno,permeab_matriz,viscosidade,propC,Nc,num_phases,deltaT) #cálculo das pressões no passo de tempo t
            phase_velocities = IMPEC.Darcy(n_phases,p,pc,D,phase_density)

            for i in range(0,n): #contador para cada bloco da malha
                for k in range(0,k): #contador para cada componente
                    Nc[k] = IMPEC.composition(propComponente[k],phase_velocities,num_phases,deltaT) #depende das propriedades do componente no passo de tempo anterior

                S = IMPEC.saturation(num_phases, p, p_old,deltaT, propC)

                propC = IMPEC.properties() #provavelmente vai ta em outra classe/outro arquivo

                estabilitity = phaseStabilityTest(propC)
                # Não tenho certeza se seria exatamente isso, acredito que só irei confirmar implementando
                while estabilitity == 0:
                    phase new = flash_calculation
                    propC = IMPEC.properties()
                    estabilitity = phaseStabilityTest(propC)

            p_old = p

            deltaT = IMPEC.times(t,tf,CFL,h,v) #duvida quanto ao valor de v
            t = t + deltaT
            #checar t<tfinal
        return P,t
## Checagem do balanço de massa

import numpy as np
from .Checa_Estabilidade import StabilityTest



class Multiphase(object):
    def __init__(self,w,Bin,R,Tc,Pc,T,P,Nc,C7,z):
        self = StabilityTest(w,Bin,R,Tc,Pc,T,P,Nc,C7,z)

    def find_fugacity_coeff(self):
        # Aqui a correlação de wilson é utilizada apenas para achar o K inicial,
        while (sum((1-fl/fv)**2))> 10e-10:

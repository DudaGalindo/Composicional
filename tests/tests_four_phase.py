import numpy as np
from stability_check import StabilityCheck
import unittest


class testes_caso_william(unittest.TestCase):

    def test1(self):
        R = 8.3144598
        z = np.array([[0.9], [0.1], [0.0]])#[n_bXn_k]
        P = np.array([7e6])
        T = np.array([311])
        Tc = np.array([190.56, 305.32, 369.83])
        Pc = np.array([4599000, 4873000, 4248000])
        Mw = np.array([16.043e-3, 30.070e-3, 44.096e-3])
        w = np.array([0.0115, 0.0995, 0.1523])
        Bin = np.array([[0.,0.,.0], [0.,0.,0.], [0.,0.,.0]])
        obj = StabilityCheck(w,Bin,R,Tc,Pc,T,P)
        obj.run(z,Mw)
        import pdb; pdb.set_trace()

#from .. import directories as direc
import numpy as np

class BrooksAndCorey:

    def __init__(self, Sorw, Sorg, Sgr, Swc, n_w, n_o, n_g, krw0, krg0, kro0):
        self.Sorw = Sorw #float(direc.data_loaded['compositional_data']['residual_saturations']['Sorw'])
        self.Sorg = Sorg #float(direc.data_loaded['compositional_data']['residual_saturations']['Sorg'])
        self.Sgr =  Sgr #float(direc.data_loaded['compositional_data']['residual_saturations']['Sgr'])
        self.Swc = Swc #float(direc.data_loaded['compositional_data']['residual_saturations']['Swc'])

        self.n_w = n_w #float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_w'])
        self.n_o = n_o #float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_o'])
        self.n_g = n_g #float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_g'])

        # End-point relative permeability data
        self.krw0 = krw0 #float(direc.data_loaded['compositional_data']['relative_permeability_data']['krw0'])
        self.kro0 = kro0 #float(direc.data_loaded['compositional_data']['relative_permeability_data']['kro0'])
        self.krg0 = krg0 #float(direc.data_loaded['compositional_data']['relative_permeability_data']['krg0'])

    def relative_permeabilities(self, saturations):
        #saturations = [So,Sg,Sw]
        self.Sor = self.Sorw * (1 - saturations[1] / (1 - self.Swr - self.Sorg)) + \
                    self.Sorg * saturations[1] / (1 - self.Swr - self.Sorg)
        kr = np.zeros(saturations.shape)
        krw = self.krw0 * ((saturations[2] - self.Swr) / (1 - self.Swr - self.Sorw - self.Sgr)) ** self.n_w
        kro = self.kro0 * ((saturations[0] - self.Sor) / (1 - self.Swr - self.Sorw - self.Sgr)) ** self.n_o
        krg = self.krg0 * ((saturations[1] - self.Sgr) / (1 - self.Swr - self.Sorw - self.Sgr)) ** self.n_g
        return kro, krg, krw

    def __call__(self, saturations):
        return self.relative_permeabilities(saturations)

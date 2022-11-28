# MANUALLY ADJUST HIPP CELL PARAMS
# -----------------------------

import matplotlib
import numpy as np

# nicer font options:
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 12})

matplotlib.use('Agg')  # hopefully this works over ssh
import matplotlib.pyplot as plt
import pylab
from random import Random  # TODO replace with numpy rand f'n.  pseudorandom number generation
from inspyred import ec  # evolutionary algorithm
from netpyne import specs, sim  # neural network design and simulation
from clamps import IClamp
from IVdata import IVdata
# from clamps_noise import ICNoise
from find_rheobase import ElectrophysiologicalPhenotype
from scipy.signal import find_peaks
from tabulate import tabulate


netparams = specs.NetParams()

hipp = netparams.importCellParams(
    label='HIPP',
    conds={"cellType": "HIPPCell", "cellModel": "HIPPCell"},
    fileName="objects/HIPP.hoc",
    cellName="HIPPCell",
    cellArgs=[1],
    importSynMechs=False
)


with open('figures/hippcell/hipp.txt', 'w') as f:
    f.write(str(hipp))


class testparam(object):
    def __init__(self,
                 cell,
                 flag):
        self.cell_dict = {"secs": cell["secs"]}
        self.flag = flag

    def curr_inj(self, current, delay=100, duration=600):
        iclamp = IClamp(self.cell_dict, delay=delay, duration=duration, T=duration + delay * 2)
        res = iclamp(current)
        return res

    def volt_inj(self):
        IV = IVdata(self.cell_dict)
        self.testclamp = IV.compute_ivdata(vlow=-70, vhigh=20, n_steps=10, delay=10, duration=5)
        return self.testclamp

    def sim_fi(self):
        ep = ElectrophysiologicalPhenotype(self.cell_dict)
        self.simfi = ep.compute_fi_curve(ilow=0, ihigh=0.4, n_steps=14, delay=0, duration=1000)
        return self.simfi

    def manual_adjust(self):
        baseline = self.sim_fi()
        baselineiv = self.volt_inj()
        
        self.cell_dict['secs']['soma']['geom']['cm'] = 1.55

        self.cell_dict['secs']['soma']['mechs']['bk']['gkbar'] = 0.0012853021551045662

        # --- SODIUM
        self.cell_dict['secs']['soma']['mechs']['ichan2']['gnatbar'] = 0.3602593904782172
        self.cell_dict['secs']['soma']['mechs']['ichan2']['vshiftma'] =  59 #44.125373392717904
        self.cell_dict['secs']['soma']['mechs']['ichan2']['vshiftmb'] = 20.0284363930537
        self.cell_dict['secs']['soma']['mechs']['ichan2']['vshiftha'] =  137.59318591477995
        self.cell_dict['secs']['soma']['mechs']['ichan2']['vshifthb'] = 11.207695065430824  # 12.5  # of interest

        # --- POTASSIUM

        self.cell_dict['secs']['soma']['mechs']['ichan2']['gkfbar'] = 0.1 #0.02472940842061752 #22   # 0.020484307 #0.016
        self.cell_dict['secs']['soma']['mechs']['ichan2']['gksbar'] = 0.002173092078840193 #0.002
        self.cell_dict['secs']['soma']['mechs']['ichan2']['vshiftnfa'] = 30 #28.193130792988576
        self.cell_dict['secs']['soma']['mechs']['ichan2']['vshiftnfb'] = 60 #47.554454544883015
        self.cell_dict['secs']['soma']['mechs']['ichan2']['vshiftnsa'] = 45 #31.539655883127935
        self.cell_dict['secs']['soma']['mechs']['ichan2']['vshiftnsb'] = 100 #86.43028799026294

        self.cell_dict['secs']['soma']['mechs']['ichan2']['gl'] = 3.053617351629548e-05

        self.cell_dict['secs']['soma']['mechs']['ka']['gkabar'] = 0.001550679735649292 #0.0008
        #self.cell_dict['secs']['soma']['mechs']['kir']['gkbar'] = 0.0
        self.cell_dict['secs']['soma']['mechs']['km']['gbar'] = 0.001
        self.cell_dict['secs']['soma']['mechs']['lca']['glcabar'] = 0.0007274383172642894 #0.00017469257304634966 #0.0015
        self.cell_dict['secs']['soma']['mechs']['nca']['gncabar'] = 0 #3.0246490314670187e-05 #0.0
        self.cell_dict['secs']['soma']['mechs']['sk']['gskbar'] = 0.0044681415100893155 #0.022809187461813662 #0.003
        #self.cell_dict['secs']['soma']['mechs']['tca']['gcatbar'] = 3.7E-05
        self.cell_dict['secs']['soma']['mechs']['ih']['ghyfbar'] = 0 #3.532543597323804e-05
        self.cell_dict['secs']['soma']['mechs']['ih']['ghysbar'] = 0 #3e-06
        
        shifted = self.sim_fi()
        shiftediv = self.volt_inj()
        self.cell_dict['secs']['soma']['mechs']['ichan2']['vshiftma'] = 43
        self.cell_dict['secs']['soma']['mechs']['ichan2']['vshiftha'] = 65
        self.cell_dict['secs']['soma']['mechs']['ichan2']['gksbar'] = 0  # 0.002
        #aptest = self.curr_inj(4)
        
        plt.plot(baseline['I'], baseline['F'], label="baseline")
        plt.plot(shifted['I'], shifted['F'], label="manual fit")
        plt.legend()
        plt.savefig("figures/hippcell/ifvsif_%s.jpeg" % self.flag)
        plt.close()

        plt.plot(baselineiv['V'], baselineiv['Na'], label='BL Na')
        plt.plot(baselineiv['V'], baselineiv['K'], label='BL K')
        plt.plot(shiftediv['V'], shiftediv['Na'], label='SH Na')
        plt.plot(shiftediv['V'], shiftediv['K'], label='SH K')
        plt.legend()
        plt.savefig("figures/hippcell/ivvsiv_%s.jpeg" % self.flag)
        plt.close()
        
        #plt.plot(aptest['t'], aptest['V'])
        #plt.savefig("figures/manual-adjust/aptest_%s.jpeg" % self.flag)


TestParam = testparam(hipp, "hipp")
TestParam.manual_adjust()

#TestParam = testparam(bc, "bc")
#TestParam.manual_adjust()
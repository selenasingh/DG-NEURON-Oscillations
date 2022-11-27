import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netpyne import sim, specs
from scipy.optimize import golden
from clamps import VClamp, IClamp

netparams = specs.NetParams()

# import granule cell info from hoc
gc = netparams.importCellParams(
    label='GC',
    conds={"cellType": "GranuleCell", "cellModel": "GranuleCell"},
    fileName="objects/GC.hoc",
    cellName="GranuleCell",
    cellArgs=[1],
    importSynMechs=False
)


class IVdata(object):
    """returns IV curve for various injected voltage vals"""

    def __init__(self, cell):
        self.cell_dict = {"secs": cell["secs"]}
        self.iv_interm = {}
        self.iv_data = None

    def step_voltage(self, voltage, delay=10, duration=5):
        """
        
        Injects voltage, returns current

        Parameters
        ----------
        voltage : FLOAT
            voltage at which membrane is clamped to [mV]
        delay : FLOAT, optional
            Time after recording starts where voltage is clamped [ms]. The default is 10.
        duration : FLOAT, optional
            Total duration of simulation [ms]. The default is 5.

        Returns
        -------
        DICT
        results of voltage clamp

        """
        vclamp = VClamp(self.cell_dict, delay=delay, duration=duration, T=duration + delay * 2)
        results = vclamp(voltage)
        return results

    def compute_ivdata(self, vlow=-70, vhigh=20, n_steps=10, delay=10, duration=5):
        self.iv_data = pd.DataFrame(np.zeros((n_steps, 3)), columns=["V", "Na", "K"])
        voltage_steps = np.linspace(vlow, vhigh, n_steps)  # creates same volt clamp steps used in experiments
        self.iv_data["V"] = voltage_steps

        for j, voltage in enumerate(voltage_steps):
            voltclamp = self.step_voltage([voltage], delay=delay, duration=duration)
            min_ina = np.min(voltclamp['i_na'][0:delay])
            max_ik = np.max(voltclamp['i_k'][0:delay])
            self.iv_data.iloc[j, 0] = voltage
            self.iv_data.iloc[j, 1] = min_ina
            self.iv_data.iloc[j, 2] = max_ik

        return self.iv_data



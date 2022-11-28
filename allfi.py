### PLOT ALL FI CURVES ###
# -----------------------
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')  # hopefully this works over ssh
import matplotlib.pyplot as plt
from netpyne import specs, sim  # neural network design and simulation
from clamps import IClamp
from find_rheobase import ElectrophysiologicalPhenotype

netparams = specs.NetParams()

mc = netparams.importCellParams(
        label='MC',
        conds={"cellType": "MossyCell", "cellModel": "MossyCell"},
        fileName="objects/MC.hoc",
        cellName="MossyCell",
        cellArgs=[1],
        importSynMechs=False
    )

gc = netparams.importCellParams(
        label='GC',
        conds={"cellType": "GranuleCell", "cellModel": "GranuleCell"},
        fileName="objects/GC.hoc",
        cellName="GranuleCell",
        cellArgs=[1],
        importSynMechs=False
    )

bc = netparams.importCellParams(
        label='BC',
        conds={"cellType": "BasketCell", "cellModel": "BasketCell"},
        fileName="objects/BC.hoc",
        cellName="BasketCell",
        cellArgs=[1],
        importSynMechs=False
    )

hipp = netparams.importCellParams(
        label='HIPP',
        conds={"cellType": "HIPPCell", "cellModel": "HIPPCell"},
        fileName="objects/HIPP.hoc",
        cellName="HIPPCell",
        cellArgs=[1],
        importSynMechs=False
    )

mc_dict = {"secs": mc["secs"]}
gc_dict = {"secs": gc["secs"]}
bc_dict = {"secs": bc["secs"]}
hipp_dict = {"secs": hipp["secs"]}



epmc = ElectrophysiologicalPhenotype(mc_dict)
mcFI = epmc.compute_fi_curve(ilow=0, ihigh=0.4, n_steps=20, delay=0, duration=1500)

epgc = ElectrophysiologicalPhenotype(gc_dict)
gcFI = epgc.compute_fi_curve(ilow=0, ihigh=0.4, n_steps=20, delay=0, duration=1500)

epbc = ElectrophysiologicalPhenotype(bc_dict)
bcFI = epbc.compute_fi_curve(ilow=0, ihigh=0.4, n_steps=20, delay=0, duration=1500)

ephipp = ElectrophysiologicalPhenotype(hipp_dict)
hippFI = ephipp.compute_fi_curve(ilow=0, ihigh=0.4, n_steps=20, delay=0, duration=1500)

plt.plot(gcFI['I'], gcFI['F'], color = 'k', linestyle = 'solid', label='gc')
plt.plot(hippFI['I'], hippFI['F'], color = 'k', linestyle = 'dashdot', label='hipp')
plt.plot(mcFI['I'], mcFI['F'], color ='k', linestyle = 'dotted', label='mc')
plt.plot(bcFI['I'], bcFI['F'], color = 'k', linestyle = 'dashed', label='bc')
plt.legend()
plt.xlabel("Current (nA)")
plt.ylabel("Frequency (Hz)")
plt.savefig('figures/mossycell/ALL_FI.jpeg')
plt.close()


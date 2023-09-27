### PLOT ALL FI CURVES ###
# -----------------------
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'  # computer modern 
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 16}) 

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
gc_dict_ma = {"secs": gc["secs"]}
bc_dict = {"secs": bc["secs"]}
hipp_dict = {"secs": hipp["secs"]}



epmc = ElectrophysiologicalPhenotype(mc_dict)
mcFI = epmc.compute_fi_curve(ilow=0, ihigh=0.6, n_steps=20, delay=0, duration=1500)

epgc = ElectrophysiologicalPhenotype(gc_dict)
gcFI = epgc.compute_fi_curve(ilow=0, ihigh=0.6, n_steps=20, delay=0, duration=1500)

#manual adjust 
gc_dict_ma['secs']['soma']['mechs']['ichan2']['gnatbar'] = (0.12 * 2.5)
gc_dict_ma['secs']['soma']['mechs']['ichan2']['gkfbar'] = (0.016 * 2.5)
gc_dict_ma['secs']['soma']['mechs']['ichan2']['gksbar'] = (0.006 * 2.5)
gc_dict_ma['secs']['soma']['mechs']['km']['gbar'] = (0.001 * 2.5)


gc_dict_ma['secs']['gcdend1_0']['mechs']['ichan2']['gnatbar'] = (0.018 * 2.5)
gc_dict_ma['secs']['gcdend1_0']['mechs']['ichan2']['gkfbar'] = (0.004 * 2.5)
gc_dict_ma['secs']['gcdend1_0']['mechs']['ichan2']['gksbar'] = (0.006 * 2.5)

gc_dict_ma['secs']['gcdend2_0']['mechs']['ichan2']['gnatbar'] = (0.018 * 2.5)
gc_dict_ma['secs']['gcdend2_0']['mechs']['ichan2']['gkfbar'] = (0.004 * 2.5)
gc_dict_ma['secs']['gcdend2_0']['mechs']['ichan2']['gksbar'] = (0.006 * 2.5)


gc_dict_ma['secs']['gcdend1_1']['mechs']['ichan2']['gnatbar'] = (0.013 * 2.5)
gc_dict_ma['secs']['gcdend1_1']['mechs']['ichan2']['gkfbar'] = (0.004 * 2.5)
gc_dict_ma['secs']['gcdend1_1']['mechs']['ichan2']['gksbar'] = (0.006 * 2.5)

gc_dict_ma['secs']['gcdend2_1']['mechs']['ichan2']['gnatbar'] = (0.013 * 2.5)
gc_dict_ma['secs']['gcdend2_1']['mechs']['ichan2']['gkfbar'] = (0.004 * 2.5)
gc_dict_ma['secs']['gcdend2_1']['mechs']['ichan2']['gksbar'] = (0.006 * 2.5)

gc_dict_ma['secs']['gcdend1_2']['mechs']['ichan2']['gnatbar'] = (0.008 * 2.5)
gc_dict_ma['secs']['gcdend1_2']['mechs']['ichan2']['gkfbar'] = (0.001 * 2.5)
gc_dict_ma['secs']['gcdend1_2']['mechs']['ichan2']['gksbar'] = (0.006 * 2.5)

gc_dict_ma['secs']['gcdend2_2']['mechs']['ichan2']['gnatbar'] = (0.008 * 2.5)
gc_dict_ma['secs']['gcdend2_2']['mechs']['ichan2']['gkfbar'] = (0.001 * 2.5)
gc_dict_ma['secs']['gcdend2_2']['mechs']['ichan2']['gksbar'] = (0.006 * 2.5)

gc_dict_ma['secs']['gcdend1_3']['mechs']['ichan2']['gnatbar'] = (0 * 2.5)
gc_dict_ma['secs']['gcdend1_3']['mechs']['ichan2']['gkfbar'] = (0.01 * 2.5)
gc_dict_ma['secs']['gcdend1_3']['mechs']['ichan2']['gksbar'] = (0.008 * 2.5)

gc_dict_ma['secs']['gcdend2_3']['mechs']['ichan2']['gnatbar'] = (0 * 2.5)
gc_dict_ma['secs']['gcdend2_3']['mechs']['ichan2']['gkfbar'] = (0.01 * 2.5)
gc_dict_ma['secs']['gcdend2_3']['mechs']['ichan2']['gksbar'] = (0.008 * 2.5)


epgc_ma = ElectrophysiologicalPhenotype(gc_dict_ma)
gcFI_ma = epgc_ma.compute_fi_curve(ilow=0, ihigh=0.6, n_steps=16, delay=0, duration=1500)

epbc = ElectrophysiologicalPhenotype(bc_dict)
bcFI = epbc.compute_fi_curve(ilow=0, ihigh=0.6, n_steps=20, delay=0, duration=1500)

ephipp = ElectrophysiologicalPhenotype(hipp_dict)
hippFI = ephipp.compute_fi_curve(ilow=0, ihigh=0.6, n_steps=20, delay=0, duration=1500)

plt.figure(figsize=(6,6))
plt.plot(gcFI['I'], gcFI['F'], color = 'k', linestyle = 'solid', label='GC')
plt.plot(gcFI_ma['I'], gcFI_ma['F'], color = '#F1586C', linestyle = 'solid', label='GC + ')
plt.plot(hippFI['I'], hippFI['F'], color = 'k', linestyle = 'dashdot', label='HIPP')
plt.plot(mcFI['I'], mcFI['F'], color ='k', linestyle = 'dotted', label='MC')
plt.plot(bcFI['I'], bcFI['F'], color = 'k', linestyle = 'dashed', label='BC')
plt.legend()
plt.xlabel("Current (nA)")
plt.ylabel("Frequency (Hz)")
plt.savefig('figures/mossycell/ALL_FI.pdf')
plt.close()

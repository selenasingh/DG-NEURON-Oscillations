# PARAMETER OPTIMIZATION TO FIX HIPP CELL FI CURVES
# --------------------
# DATA USED FOR FITTING:
# Degro CE, Bolduan F, Vida I, Booker SA. Interneuron diversity in the rat dentate gyrus: An unbiased in vitro 
# classification. Hippocampus. 2022;32(4):310–331. doi:10.1002/hipo.23408
# --------------------
# NOTE: Parameters for SOMA ONLY are optimized here, to simplify re-implementation in NEURON
#       Fit isn't great, but better than baseline. 


import matplotlib
import numpy as np

# nicer font options:
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 12})

matplotlib.use('Agg')  # hopefully this works over ssh
import matplotlib.pyplot as plt
import pandas as pd
from random import Random
from inspyred import ec  # evolutionary algorithm
from netpyne import specs, sim  # neural network design and simulation
from clamps import IClamp
from find_rheobase import ElectrophysiologicalPhenotype

netparams = specs.NetParams()

hipp = netparams.importCellParams(
    label='HIPP',
    conds={"cellType": "HIPPCell", "cellModel": "HIPPCell"},
    fileName="objects/HIPP.hoc",
    cellName="HIPPCell",
    cellArgs=[1],
    importSynMechs=False
)

# parameters to be optimized
free_params = {
    'bk': ['gkbar'],  # big conductance, calcium-activated potassium channel
    'ichan2': ['gnatbar', 'vshiftma', 'vshiftmb', 'vshiftha', 'vshifthb', 'vshiftnfa', 'vshiftnfb', 'vshiftnsa',
               'vshiftnsb',
               'gkfbar', 'gksbar', 'gl'],  # sodium, potassium parameters
    'ka': ['gkabar'],  # A-type (fast inactivating) Kv channel
    'lca': ['glcabar'],  # l-type calcium
    'nca': ['gncabar'],  # n-type calcium
    'sk': ['gskbar'],  # small conductance potassium channel
    'ih': ['ghyfbar', 'ghysbar']  # HCN channel
}


class testparam(object):
    def __init__(self,
                 cell,
                 free_params,
                 pop_size=10,
                 max_evaluations=60,
                 num_selected=10,
                 mutation_rate=0.03,
                 ):
        self.cell_dict = {"secs": cell["secs"]}
        self.free_params = free_params
        self.initialParams = []
        self.minParamValues = []
        self.maxParamValues = []
        self.num_inputs = len(sum(self.free_params.values(), []))
        self.free_params = free_params
        self.pop_size = pop_size
        self.max_evaluations = max_evaluations
        self.num_selected = num_selected
        self.mutation_rate = mutation_rate
        self.num_elites = 1
        self.flag = str('hipp')
        self.n_simcells = 1  # number of simulated cells

        # self.manual_adjust()
        self.plot_results()

    def retrieve_baseline_params(self):
        """ Saves baseline parameters from cell_dict

        Returns:
            'list'. List of baseline parameter values.
        """
        self.baseline = []
        for key in self.free_params.keys():
            for val in self.free_params[key]:
                self.baseline.append(self.cell_dict['secs']['soma']['mechs'][key][val])
        self.num_inputs = len(self.baseline)
        return self.baseline

    def curr_inj(self, current, delay=0, duration=1000):
        iclamp = IClamp(self.cell_dict, delay=delay, duration=duration, T=duration + delay * 2)
        res = iclamp(current)
        return res

    def sim_fi(self):
        ep = ElectrophysiologicalPhenotype(self.cell_dict)
        self.simfi = ep.compute_fi_curve(ilow=0, ihigh=0.25, n_steps=6, delay=0, duration=1000)
        return self.simfi

    def data_fi(self):
        x = [0, 0.05, 0.1, 0.15, 0.2, 0.25]
        y = [0, 7, 19, 40, 51, 61]
        datafi = [x, y]
        self.datafi = np.array(datafi)
        return self.datafi

    def generate_netparams(self, random, args):
        """
        Initialize set of random initial parameter values selected from uniform distribution within min-max bounds.

        Returns
            'list'. initialParams
        """
        self.initialParams = [random.uniform(self.minParamValues[i], self.maxParamValues[i]) for i in
                              range(self.num_inputs)]
        return self.initialParams

    def evaluate_netparams(self, candidates, args):
        """
        Fitness function that evaluates the fitness of candidate parameters by quantifying the difference between
        simulated FI and IV curves to the FI and IV curves from data using mean squared error.

        Returns
            'list'. Fitness values for sets of candidates
        """
        self.fitnessCandidates = []

        for cand in candidates:
            i = 0
            for k in free_params.keys():
                for v in free_params[k]:
                    self.cell_dict['secs']['soma']['mechs'][k][v] = cand[i]
                    i += 1
            FI_data = self.data_fi()
            FI_sim = self.sim_fi().to_numpy()

            ficurves = np.sum([((x1 - x2) ** 2) for (x1, x2) in zip(FI_data[1, :], FI_sim[:, 1])]) / len(FI_data[1, :])

            fitness = ficurves

            self.fitnessCandidates.append(fitness)

        return self.fitnessCandidates

    def find_bestcandidate(self):
        """
        Sets up custom evolutionary computation and returns list of optimized parameters.

        Components of EC
            gc_ec : instantiate evolutionary computation with random.Random object
            selector: method used to select best candidate based on fitness value
            variator: method used to determine how mutations (variations) are made to each generation of params
            replacer: method used to determine if/how individuals are replaced in pool of candidates after selection
            terminator: method used to specify how to terminate evolutionary algorithm
            observer: method that allows for supervision of evolutionary computation across all evaluations
            evolve: pulls together all components of custom algorithm, iteratively calls evaluate_netparams, returns
                    parameters that minimize fitness function.

        Returns
            'list'. bestCand (list of optimized parameters)
        """

        # TODO: Potentially write custom variator function to be compatible with np.random.RandomState
        # rand = np.random.RandomState(self.setseed)

        rand = Random()
        rand.seed(self.setseed)  # will take cell # as seed (n_simcells = 1, seed = 0. n_simcells = 2, seed = 1, etc).

        # SET UP MIN/MAX BOUNDS FOR PARAMETERS ------------------
        # TODO: find cleaner way of dealing with these lists, allow for easier modification

        #self.minParamValues = [0.2 * param for param in self.baseline]  # 0.5 best for IF, Na
        #self.maxParamValues = [3.0 * param for param in self.baseline]  # 3.0 best for IF, Na

        
        self.minParamValues = [(0.0006 * 0.2), (0.2 * 1.0), (43*0.9), (15 * 0.2), (100*0.2), (12.5 * 0.2),
                               (18.0 * 0.2), (43.0 * 0.2), (30.0 * 0.2), (55.0 * 0.2), (0.01 * 0.2), (0.01 * 0.9), (3.6e-05 * 0.2),
                               (0.0008* 0.2), (0.0015 * 0.2), (8e-05 * 0.2), (0.003 * 0.2),
                               (1.5e-05 * 0.2), (1.5e-05 * 0.2)
                               ]

        self.maxParamValues = [(0.0006 * 2.5), (0.2 * 2.5), (43*1.2), (15 * 2.5), (100*2.5), (12.5 * 2.5),
                               (18.0 * 2.5), (43.0 * 2.5), (30.0 * 2.5), (55.0 * 2.5), (0.01 * 2.5), (0.01 * 10.0), (3.6e-05 * 2.5),
                               (0.0008 * 2.5), (0.0015 * 2.5), (8e-05 * 2.5), (0.003 * 2.5),
                               (1.5e-05 * 2.5), (1.5e-05 * 2.5)
                               ]
        
        # SET UP EVOLUTIONARY COMPUTATION ----------------------
        self.gc_ec = ec.EvolutionaryComputation(rand)
        self.gc_ec.selector = ec.selectors.truncation_selection  # purely deterministic
        self.gc_ec.variator = [ec.variators.uniform_crossover, ec.variators.gaussian_mutation]
        self.gc_ec.replacer = ec.replacers.generational_replacement
        self.gc_ec.terminator = ec.terminators.evaluation_termination  # terminates after max number of evals is met
        self.gc_ec.observer = ec.observers.plot_observer  # save to file, use observers.file_observer

        self.final_pop = self.gc_ec.evolve(generator=self.generate_netparams,  # f'n for initializing params
                                           evaluator=self.evaluate_netparams,  # f'n for evaluating fitness values
                                           pop_size=self.pop_size,  # number of parameter sets per evaluation
                                           maximize=False,  # best fitness corresponds to minimum value
                                           bounder=ec.Bounder(  # set min/max param bounds
                                               self.minParamValues,
                                               self.maxParamValues
                                           ),
                                           max_evaluations=self.max_evaluations,
                                           num_selected=self.num_selected,
                                           mutation_rate=self.mutation_rate,
                                           num_inputs=self.num_inputs,
                                           num_elites=self.num_elites
                                           )

        self.final_pop.sort(reverse=True)  # sort final population so best fitness is first in list
        self.bestCand = self.final_pop[0].candidate  # bestCand <-- individual @ start of list

        plt.savefig('figures/hippcell/observer_%s.pdf' % self.flag)  # save fitness vs. iterations graph
        plt.close()
        return self.bestCand

    def build_optimizedcell(self):
        """ Replaces baseline parameters with parameters from best candidate, then uses current injection experiment
            to build 'optimized' cell.

        Returns
            'dict'. Results of current clamp from optimized cell.
        """
        j = 0
        for key in self.free_params.keys():
            for val in self.free_params[key]:
                self.cell_dict['secs']['soma']['mechs'][key][val] = self.bestCand[j]
                j = j + 1
        finalclamp = self.curr_inj(0.33)
        # self.ep_opt = ElectrophysiologicalPhenotype(self.cell_dict)
        return finalclamp

    def store_curves(self):
        """ Generates set of n optimized neurons (n_simcells), stores baseline and optimized parameters,
            IF and IV curves.

        Returns
            'tuple' of two DataFrames, (sim_fi_store, sim_iv_store)
        """
        # initialize empty DataFrames, populate with baseline parameters
        baselineparams = self.retrieve_baseline_params()
        self.param_store = pd.DataFrame({"param": sum(free_params.values(), []), "baseline": baselineparams})
        self.sim_fi_store = pd.DataFrame([])

        # generate set of n_simcells, populate DataFrames above with FI, IV, params
        for cell_n in range(0, self.n_simcells):
            self.setseed = cell_n  # set new seed for evol'n computation
            newparams = self.find_bestcandidate()  # find optimized parameters
            newparamdf = pd.DataFrame({"Cell_%s" % cell_n: newparams})  # store those params with a label
            self.param_store = pd.concat([self.param_store, newparamdf], axis=1)  # append params to DF
            self.build_optimizedcell()  # build the optimized cell
            newcellfi = self.sim_fi()  # generate simulated FI curve
            self.sim_fi_store = pd.concat([newcellfi, self.sim_fi_store])  # append FI curve to DF

        # save dataframes to .csv
        self.sim_fi_store.to_csv('data/parameters/simFIs_%s.csv' % self.flag)
        self.param_store.to_csv('data/parameters/parameters_%s.csv' % self.flag)

        return self.sim_fi_store 

    def compute_avg_curves(self):
        """ Computes average simulated FI and IV curves and SEM

        Returns
            'tuple' of two DataFrames, (avg_FI, avg_IV)
        """
        sim_fi_store = self.store_curves()

        # average simulated FI curve:
        avgfi = sim_fi_store.groupby(['I']).agg({'F': ['mean']}).values
        semfi = sim_fi_store.groupby(['I']).agg({'F': ['std']}).values / np.sqrt(self.n_simcells)
        self.avg_FI = np.c_[np.linspace(0, 0.25, 6), avgfi, semfi]

        return self.avg_FI 

    def plot_results(self):
        """ Plots average simulated IV and FI curves from optimized neurons against avg curves from data. Saves
        figure to 'figures/hippcell'. Automatically called when optimizeparams is instantiated.
        """

        # Generate and collect all data for plotting
        baselinecellfi = self.sim_fi().to_numpy()
        exp_fi = self.data_fi()
        avg_fi = self.compute_avg_curves()

        fig1, (ax1) = plt.subplots(1, 1)
        
        # FI curves
        ax1.plot(baselinecellfi[:, 0], baselinecellfi[:, 1], color='0.7', linestyle='dashed', label='Baseline')
        ax1.plot(exp_fi[0, :], exp_fi[1, :], color='0.5', label='Data')
        ax1.errorbar(avg_fi[:, 0], avg_fi[:, 1], yerr=avg_fi[:, 2], color='0.0', label='Optimized')
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax1.set_xlabel("Current (nA)")
        ax1.set_ylabel("Frequency (Hz)")

        fig1.tight_layout()
        fig1.savefig('figures/hippcell/optimizationresults_%s.pdf' % self.flag, bbox_inches="tight")


TestParam = testparam(hipp, free_params)

import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import gaussian_kde as kde
#from progressbar import Percentage,Bar,RotatingMarker,ETA,ProgressBar
import re
from scipy.integrate import quad
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from scipy.stats import gaussian_kde


import matplotlib.patheffects as PathEffects
import matplotlib.gridspec as gridspec
import glob
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times','Palatino']})
rc('text', usetex=True)

from dmdd import *
from experiment_parameters import *
from multinest_parameters import *

mpl.rcParams['xtick.major.size']=8
mpl.rcParams['ytick.major.size']=8
mpl.rcParams['xtick.labelsize']=18
mpl.rcParams['ytick.labelsize']=18

try:
    MAIN_PATH = os.environ['DMDD_AM_MAIN_PATH']
except KeyError:
    logging.warning('DMDD_AM_MAIN_PATH environment variable not defined, defaulting to current working directory.')
    MAIN_PATH = os.getcwd()

SIM_PATH = MAIN_PATH + '/simulations_uv/'
CHAINS_PATH = MAIN_PATH + '/chains_uv/'
RESULTS_PATH = MAIN_PATH + '/results_uv'

if not os.path.exists(SIM_PATH):
    os.mkdir(SIM_PATH)
if not os.path.exists(CHAINS_PATH):
    os.mkdir(CHAINS_PATH)
if not os.path.exists(RESULTS_PATH):
    os.mkdir(RESULTS_PATH)

###############################################

SI_Higgs = UV_Model('SI_Higgs', ['mass', 'sigma_si'],
                    fixed_params={'fnfp_si': 1}, time_info=True,GF=True)
anapole = UV_Model('Anapole', ['mass','sigma_anapole'], time_info=True,GF=True)
ALLMODELS = [SI_Higgs, anapole]

SI_Higgs_noT = UV_Model('SI_Higgs', ['mass', 'sigma_si'],
                    fixed_params={'fnfp_si': 1}, time_info=False,GF=True)
anapole_noT = UV_Model('Anapole', ['mass','sigma_anapole'], time_info=False,GF=True)
ALLMODELS = [SI_Higgs, anapole]


################################################################################################

experiment = 'Xe'
xe = dmdd.Experiment(experiment,Target[experiment],
                     Qmin[experiment], Qmax[experiment],
                     Exposure[experiment], eff.efficiency_Xe,
                     Start[experiment], End[experiment])

experiment = 'Xe10x'
xe10x = dmdd.Experiment(experiment,Target[experiment],
                       Qmin[experiment], Qmax[experiment],
                       Exposure[experiment], dmdd.eff.efficiency_Xe,
                       Start[experiment], End[experiment])


ALL_EXPERIMENTS = {}
ALL_EXPERIMENTS['Xe'] = xe
ALL_EXPERIMENTS['Xe10x'] = xe10x
#ALL_EXPERIMENTS['Ge+Xe'] = [ge, xe]

################################################################################################
sigma_names = {}
fnfp_names = {}
fnfp_vals = {}

for m in ALLMODELS:        
    sigma_names[m.name] = m.param_names[1]
    
    if len(m.fixed_params) > 0:
        fnfp_names[m.name] = m.fixed_params.keys()[0]
        fnfp_vals[m.name] = m.fixed_params.values()[0]
    else:
        fnfp_names[m.name] = None
        fnfp_vals[m.name] = None


mmin = 1
mmax = 1000

MASSES=np.logspace(np.log10(mmin),np.log10(mmax),100)


################################################################################################

def line_plots_1vsall(sigma_lim_file,
                      nsim=1, startsim=1, masses=[20.],
                      experiment_names=['Xe10x'],
                      simmodels=[SI_Higgs], models=[SI_Higgs_noT,anapole_noT],
                      filelabel='withT_withGF', allverbose=False, verbose=False,
                      results_root='/data/verag/dmdd-am/results_uv/',
                      saveplots=True, alpha=0.3, xoffset=0.1,
                      fs=20, fs2=18,
                      analyze_energy=True):
  
    if len(filelabel)>0:
        filelabel = '_' + filelabel
   
    sigma_limvals = pickle.load(open(sigma_lim_file,'rb'))
    experiment_labels = []#experiment_names
    experiments = []
    for i,exn in enumerate(experiment_names):
        experiments.append(ALL_EXPERIMENTS[exn])
        experiment_labels.append(Experiment_LaTeX[exn])


        
    for mass in masses:
        print '{}GeV:'.format(mass)

        
        for m in simmodels:
            print ''
            print ''
            print m.name
            plt.figure()
            ax=plt.gca()
            ticks = np.arange(len(experiments)) + 0.5
            
            for i,experiment in enumerate(experiments):
                try:
                    if len(experiment) > 1:
                        experimentlist = experiment
                except:
                    experimentlist = [experiment]
                ys = np.zeros(nsim)
                
                for ni,n in enumerate(np.arange(startsim, nsim + startsim)):
                    renorm = -np.inf
                    ev_sum = 0.
                    ln_ev = np.zeros(len(models))
                    prob = np.zeros(len(models))
                    
                    for j,fitm in enumerate(models):
                        pardic = {'mass': mass, sigma_names[m.name]: sigma_limvals[mass][m.name]}
                        mnrun = MultinestRun('sim{}'.format(n), experimentlist, m,
                                             pardic,fitm, prior_ranges=prior_ranges,
                                            force_sim=False,n_live_points=n_live_points,
                                            silent=True,empty_run=True,
                                            analyze_energy=analyze_energy)
                        if allverbose:
                            print ''
                            print mnrun.foldername
                        if verbose and not allverbose and fitm.name == m.name:
                            print ''
                            print '({})'.format(mnrun.foldername)
                            
                        ln_ev[j] = mnrun.get_evidence()
                        renorm = np.amax(np.array([renorm,ln_ev[j]]))
                        
                    for j,fitm in enumerate(models):
                        ev_sum += np.exp(-renorm + ln_ev[j])
                        
                    for j,fitm in enumerate(models):
                        prob[j] = np.exp(-renorm + ln_ev[j]) / ev_sum
                        if fitm.name != m.name:
                            if verbose or allverbose:
                                print '{} {}'.format(fitm.name, prob[j])
                        else:
                            ys[ni] = prob[j]
                            if verbose or allverbose:
                                print '{} {} <--sim'.format(fitm.name, prob[j])
                            if verbose and prob[j]<0.1:
                                print 'FLAG! truth unlikely?   <-----------|'
                                print ''
                                

                    if verbose:
                        print ''
                        
                    if ev_sum == np.inf:
                        print 'yikes. evidence is infinity.'
                    
                    
                ymedian = np.median(ys)
                ylo, yhi = np.percentile(ys, 25), np.percentile(ys, 75)
                for y in ys:
                    plt.plot( [ i + xoffset, i + 1 - xoffset ], [ y, y ], lw=1, 
                              alpha=alpha, color=Colors[experiment_names[i]])
                plt.plot( [ i + 0.5, i + 0.5 ], [ ylo, yhi ], color='k', lw=2)
                plt.plot( i + 0.5, ymedian, 'x', ms = 10, color='k', mew=2 )
                success = float(np.sum( ys > 0.9 )) / len( ys ) * 100.

                if success > 0:
                    plt.annotate( '{:.0f}\\%%'.format(success), xy = ( i + 0.5, 1.05 ), 
                                  fontsize=fs2, va='center', ha='center', color='k' )
                    
            ax.set_xticks( ticks )
            ax.set_xticklabels( experiment_labels )
            ax.set_title('True model: {} (mass: {:.0f} GeV)'.format(MODELNAME_TEX[m.name], mass), fontsize=fs)
            pl.ylim( ( -0.04, 1.19 ) )
            pl.axhline(0.9, ls='--', color='k', lw=2)
          
            pl.annotate('Simulations preferring true model (with \\textgreater 0.9 probability):', 
                        xy=( 0.02, 0.96 ), xycoords='axes fraction', fontsize=fs2, va='top', color='k')
            pl.ylabel('Probability of true model', fontsize=fs)

            figname = results_root + 'lineplot_{:.0f}GeV_{}_{}sims{}.pdf'.format(mass, m.name, nsim,filelabel)
            if saveplots:
                pl.savefig(figname)
            
########################

color_groups = [('Aqua', 'Blue'),
                ('Orange', 'Red'),
                ('Magenta', 'Pink')]
##########################################
def get_chains(sim_name, experiments, sim_model,
            param_values, fit_model, 
            prior='logflat',
            sim_root=SIM_PATH, chains_root=CHAINS_PATH,
            n_live_points=2000, basename='1-',
            params_to_get=None):

    foldername = sim_name
    experiment_names = [ex.name for ex in experiments]
    experiment_names.sort()
    for experiment in experiment_names:
        foldername += '_{}'.format(experiment)

    all_params = dict(param_values)
    for k in sim_model.fixed_params:
        all_params[k] = sim_model.fixed_params[k]

    all_param_names = sim_model.param_names + sim_model.fixed_params.keys()
    inds = np.argsort(all_param_names)        
    sorted_param_names = np.array(all_param_names)[inds]
    sorted_param_values = [all_params[par] for par in sorted_param_names]
    for parname, parval in zip(sorted_param_names, sorted_param_values):
        foldername += '_{}_{:.2f}'.format(parname, parval)      
    foldername += '_fitfor'

    inds = np.argsort(fit_model.param_names)
    sorted_fitparam_names = np.array(fit_model.param_names)[inds]
    for par in sorted_fitparam_names:
        foldername += '_{}'.format(par)

    if len(fit_model.fixed_params) > 0:
        foldername += '_fixed'

        keys = fit_model.fixed_params.keys()
        inds = np.argsort(keys)
        sorted_fixedparam_names = np.array(keys)[inds]
        sorted_fixedparam_values = [fit_model.fixed_params[par] for par in sorted_fixedparam_names]
        for parname, parval in zip(sorted_fixedparam_names, sorted_fixedparam_values):
            foldername += '_{}_{:.2f}'.format(parname, parval)   

    foldername += '_{}_nlive{}'.format(prior, n_live_points)
    
    chainspath = '{}/{}/'.format(chains_root, foldername)

    chainsfile = chainspath + '/' + basename + 'post_equal_weights.dat'

    data = np.loadtxt(chainsfile)

    if params_to_get is None:
        params_to_get = fit_model.param_names
    elif type(params_to_get) != list:
        params_to_get = [params_to_get]

    n_desired = len(params_to_get)
    chains = np.zeros((len(data),len(params_to_get)))
    for i,p in enumerate(params_to_get):
        ind = fit_model.param_names.index(p)
        chains[:,i] = data[:,ind]

    return chains.squeeze()

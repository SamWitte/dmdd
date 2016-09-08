import numpy as np
import matplotlib
matplotlib.use('Agg')
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
from scipy.stats import gaussian_kde, norm
import copy


import matplotlib.patheffects as PathEffects
import matplotlib.gridspec as gridspec
import glob
import dmdd
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

SIM_PATH = MAIN_PATH + '/simulations_uv'
CHAINS_PATH = MAIN_PATH + '/chains_uv'
RESULTS_PATH = MAIN_PATH + '/results_uv'

if not os.path.exists(SIM_PATH):
    os.mkdir(SIM_PATH)
if not os.path.exists(CHAINS_PATH):
    os.mkdir(CHAINS_PATH)
if not os.path.exists(RESULTS_PATH):
    os.mkdir(RESULTS_PATH)

###############################################

SI_Higgs = UV_Model('SI_Higgs', ['mass', 'sigma_si'], fixed_params={'fnfp_si': 1})
millicharge = UV_Model('Millicharge', ['mass', 'sigma_si_massless'], fixed_params={'fnfp_si_massless': 0})

SD_flavoruniversal = UV_Model('SD_fu', ['mass','sigma_sd'], fixed_params={'fnfp_sd': -1.1})
SD_Zmediated = UV_Model('SD_Z-mediated', ['mass','sigma_sd'], fixed_params={'fnfp_sd': -13.2})
SD_Moira = UV_Model('SI_Moira', ['mass','sigma_sd'], fixed_params={'fnfp_sd': 0.})

anapole = UV_Model('Anapole', ['mass','sigma_anapole'])
magdip_heavy = UV_Model('Mag.dip.heavy', ['mass','sigma_magdip'])
magdip_0 = UV_Model('Mag.dip.light', ['mass','sigma_magdip_massless'])
elecdip_heavy = UV_Model('Elec.dip.heavy', ['mass','sigma_elecdip'])
elecdip_0 = UV_Model('Elec.dip.light', ['mass','sigma_elecdip_massless'])

f1 = dmdd.UV_Model('f1', ['mass','sigma_f1'], fixed_params={'fnfp_f1': 1.})
f2_Higgs = dmdd.UV_Model('f2_Higgs', ['mass','sigma_f2'], fixed_params={'fnfp_f2': -0.05})
f2_flavoruniversal = dmdd.UV_Model('f2_flavor-universal', ['mass','sigma_f2'], fixed_params={'fnfp_f2': 1.})
f3_Higgs = dmdd.UV_Model('f3_Higgs', ['mass','sigma_f3'], fixed_params={'fnfp_f3': -0.05})
f3_flavoruniversal = dmdd.UV_Model('f3_flavor-universal', ['mass','sigma_f3'], fixed_params={'fnfp_f3': 1.})
LS = dmdd.UV_Model('LS', ['mass','sigma_LS'], fixed_params={'fnfp_LS': 0.})


MODELS_UV1 = [SI_Higgs, millicharge, SD_flavoruniversal, anapole, magdip_heavy, magdip_0, elecdip_heavy, elecdip_0]
MODELS_UV2 = [magdip_heavy, elecdip_heavy, f1, f2_Higgs, f2_flavoruniversal, f3_Higgs, f3_flavoruniversal, LS]

MODEL_CLASSES = [[SI_Higgs, SD_flavoruniversal, anapole], [magdip_heavy, elecdip_heavy], [magdip_0, elecdip_0, millicharge]]
MODEL_CLASSES_NAMES = [[SI_Higgs.name, SD_flavoruniversal.name, anapole.name], [magdip_heavy.name, elecdip_heavy.name], [magdip_0.name, elecdip_0.name, millicharge.name]]

################################################################################################

experiment = 'Xe'
xe = dmdd.Experiment(experiment,Target[experiment],
                     Qmin[experiment], Qmax[experiment],
                     Exposure[experiment], dmdd.eff.efficiency_Xe,
                     Start[experiment], End[experiment])
experiment = 'Xelo'
xelo = dmdd.Experiment(experiment,Target[experiment],
                     Qmin[experiment], Qmax[experiment],
                     Exposure[experiment], dmdd.eff.efficiency_Xe,
                       Start[experiment], End[experiment])
experiment = 'Xehi'
xehi = dmdd.Experiment(experiment,Target[experiment],
                     Qmin[experiment], Qmax[experiment],
                     Exposure[experiment], dmdd.eff.efficiency_Xe,
                       Start[experiment], End[experiment])
experiment = 'Xewide'
xewide = dmdd.Experiment(experiment,Target[experiment],
                     Qmin[experiment], Qmax[experiment],
                     Exposure[experiment], dmdd.eff.efficiency_Xe,
                         Start[experiment], End[experiment])
experiment = 'XeG3'
xeG3 = dmdd.Experiment(experiment,Target[experiment],
                     Qmin[experiment], Qmax[experiment],
                     Exposure[experiment], dmdd.eff.efficiency_Xe,
                       Start[experiment], End[experiment])
experiment = 'Ilo'
iodlo = dmdd.Experiment(experiment,Target[experiment],
                     Qmin[experiment], Qmax[experiment],
                     Exposure[experiment], dmdd.eff.efficiency_Xe,
                        Start[experiment], End[experiment])
experiment = 'Ge'
ge = dmdd.Experiment(experiment,Target[experiment],
                     Qmin[experiment], Qmax[experiment],
                     Exposure[experiment], dmdd.eff.efficiency_Xe,
                     Start[experiment], End[experiment])

experiment = 'F'
flu = dmdd.Experiment(experiment,Target[experiment],
                     Qmin[experiment], Qmax[experiment],
                     Exposure[experiment], dmdd.eff.efficiency_Xe,
                      Start[experiment], End[experiment],
                      energy_resolution=False)

experiment = 'Na'
na = dmdd.Experiment(experiment,Target[experiment],
                     Qmin[experiment], Qmax[experiment],
                     Exposure[experiment], dmdd.eff.efficiency_Xe,
                     Start[experiment], End[experiment])

experiment = 'I'
iod = dmdd.Experiment(experiment,Target[experiment],
                     Qmin[experiment], Qmax[experiment],
                     Exposure[experiment], dmdd.eff.efficiency_Xe,
                      Start[experiment], End[experiment])

experiment = 'Ar'
ar = dmdd.Experiment(experiment,Target[experiment],
                     Qmin[experiment], Qmax[experiment],
                     Exposure[experiment], dmdd.eff.efficiency_Xe,
                     Start[experiment], End[experiment])

experiment = 'F+'
flup = dmdd.Experiment(experiment,Target[experiment],
                     Qmin[experiment], Qmax[experiment],
                     Exposure[experiment], dmdd.eff.efficiency_Xe,
                     Start[experiment], End[experiment],
                     energy_resolution=False)

experiment = 'I+'
iodp = dmdd.Experiment(experiment,Target[experiment],
                     Qmin[experiment], Qmax[experiment],
                     Exposure[experiment], dmdd.eff.efficiency_Xe,
                     Start[experiment], End[experiment],
                     energy_resolution=False)

experiment = 'XeDouble'
xeplus = dmdd.Experiment(experiment,Target[experiment],
                       Qmin[experiment], Qmax[experiment],
                       Exposure[experiment], dmdd.eff.efficiency_Xe,
                         Start[experiment], End[experiment])

experiment = 'XeTriple'
xetrips = dmdd.Experiment(experiment,Target[experiment],
                       Qmin[experiment], Qmax[experiment],
                       Exposure[experiment], dmdd.eff.efficiency_Xe,
                          Start[experiment], End[experiment])

experiment = 'Xe10x'
xe10x = dmdd.Experiment(experiment,Target[experiment],
                       Qmin[experiment], Qmax[experiment],
                       Exposure[experiment], dmdd.eff.efficiency_Xe,
                       Start[experiment], End[experiment])

experiment = 'XeT'
xeT = dmdd.Experiment(experiment,Target[experiment],
                      Qmin[experiment], Qmax[experiment],
                      Exposure[experiment], dmdd.eff.efficiency_Xe,
                      Start[experiment], End[experiment])


ALL_EXPERIMENTS = {}
ALL_EXPERIMENTS['Xe'] = xe
ALL_EXPERIMENTS['Xelo'] = xelo
ALL_EXPERIMENTS['Xehi'] = xehi
ALL_EXPERIMENTS['Xewide'] = xewide
ALL_EXPERIMENTS['XeG3'] = xeG3
ALL_EXPERIMENTS['XeG3+Ilo+F'] = [xeG3, iodlo, flu]
ALL_EXPERIMENTS['Ilo'] = iodlo
ALL_EXPERIMENTS['Ge'] = ge
ALL_EXPERIMENTS['I'] = iod
ALL_EXPERIMENTS['Na'] = na
ALL_EXPERIMENTS['Ar'] = ar
ALL_EXPERIMENTS['F'] = flu
ALL_EXPERIMENTS['Ge+Xe'] = [ge, xe]
ALL_EXPERIMENTS['Ge+Xe+Na+I'] = [ge, xe, na, iod]
ALL_EXPERIMENTS['Ge+Xe+I'] = [ge, xe, iod]
ALL_EXPERIMENTS['Ge+Xe+F'] = [ge, xe, flu]
ALL_EXPERIMENTS['Ge+Xe+I+F'] = [ge, xe, iod, flu]
ALL_EXPERIMENTS['Na+I'] = [na, iod]
ALL_EXPERIMENTS['XeG3+I++F+'] = [xeG3, iodp, flup]
ALL_EXPERIMENTS['F+'] = flup
ALL_EXPERIMENTS['I+'] = iodp
ALL_EXPERIMENTS['F+Xe'] = [flu, xe]
ALL_EXPERIMENTS['XeDouble'] = xeplus
ALL_EXPERIMENTS['XeTriple'] = xetrips
ALL_EXPERIMENTS['Xe10x'] = xe10x
ALL_EXPERIMENTS['XeT'] = xeT


pandax=dmdd.Experiment('PandaX','xenon', 4., 30., 188., dmdd.eff.efficiency_Xe, 0., 1., energy_resolution=True)
lux=dmdd.Experiment('LUX','xenon', 4., 30., 30.7, dmdd.eff.efficiency_Xe, 0., 1., energy_resolution=True)
cdmslite=dmdd.Experiment('CDMSlite', 'germanium', 0.840, 6, 0.0164, dmdd.eff.efficiency_Xe,
                         0.,1.,energy_resolution=True)
supercdms=dmdd.Experiment('SuperCDMS','germanium', 1.6, 12, 1.581, dmdd.eff.efficiency_Xe,
                          0.,1.,energy_resolution=True)
kimsiod=dmdd.Experiment('KIMSiod','iodine',28.57,85.71, 32.8, dmdd.eff.efficiency_KIMS,0.,1.)#Q=0.07
kimsiodLO=dmdd.Experiment('KIMSiodLO','iodine',28.57,42.86, 32.8, dmdd.eff.efficiency_KIMS,0.,1.)
kimsiodHI=dmdd.Experiment('KIMSiodHI','iodine',42.86,300, 32.8, dmdd.eff.efficiency_KIMS,0.,1.)
picoA=dmdd.Experiment('PICO2LA','fluorine', 3.2, 400, 0.205, dmdd.eff.efficiency_Xe,0.,1.)
picoB=dmdd.Experiment('PICO2LB','fluorine', 4.4, 400, 0.046, dmdd.eff.efficiency_Xe,0.,1.)
picoC=dmdd.Experiment('PICO2LC','fluorine', 6.1, 400, 0.225, dmdd.eff.efficiency_Xe,0.,1.)
picoD=dmdd.Experiment('PICO2LD','fluorine', 8.1, 400, 0.104, dmdd.eff.efficiency_Xe,0.,1.)

current_experiments = [lux,supercdms,cdmslite,kimsiodLO,kimsiodHI,picoA,picoB,picoC,picoD]

################################################################################################
sigma_names = {}
fnfp_names = {}
fnfp_vals = {}

for m in MODELS_UV1 + MODELS_UV2:
    sigma_names[m.name] = m.param_names[1]

    if len(m.fixed_params) > 0:
        fnfp_names[m.name] = m.fixed_params.keys()[0]
        fnfp_vals[m.name] = m.fixed_params.values()[0]
    else:
        fnfp_names[m.name] = None
        fnfp_vals[m.name] = None



Nbg = {
    lux.name: 3.6,
    pandax.name: 13.95,
    supercdms.name: 4,
    cdmslite.name: 4,
    kimsiod.name: 1,
    kimsiodLO.name: 22,#500,
    kimsiodHI.name: 63,#4000,
    picoA.name: 9,
    picoB.name: 1,
    picoC.name: 3,
    picoD.name: 1,
}

mx_guess = {
    pandax.name: 1.,
    lux.name: 5.,
    supercdms.name: 5.,
    cdmslite.name: 1.,
    kimsiod.name: 1.,
    kimsiodLO.name: 1.,
    kimsiodHI.name: 1.,
    picoA.name: 1.,
    picoB.name: 1.,
    picoC.name: 1.,
    picoD.name: 1.,
}

sigma_guess = {}
for m in MODELS_UV1 + MODELS_UV2:
    sigma_guess[m.name] = 1e10

sigma_lims = {}
sigma_pandax={}
sigma_lux={}
sigma_cdmslite={}
sigma_supercdms={}
sigma_kimsiodLO={}
sigma_kimsiodHI={}
sigma_picoA={}
sigma_picoB={}
sigma_picoC={}
sigma_picoD={}

mmin = 1
mmax = 1000


MASSES=np.logspace(np.log10(mmin),np.log10(mmax),100)

for m in MODELS_UV1 + MODELS_UV2:
    sigma_lims[m.name] = np.zeros(len(MASSES))
for i,mass in enumerate(MASSES):

    sigma_lux[mass]={}
    sigma_pandax[mass] = {}
    sigma_cdmslite[mass]={}
    sigma_supercdms[mass]={}
    sigma_kimsiodLO[mass]={}
    sigma_kimsiodHI[mass]={}
    sigma_picoA[mass]={}
    sigma_picoB[mass]={}
    sigma_picoC[mass]={}
    sigma_picoD[mass]={}

    for m in MODELS_UV1 + MODELS_UV2:
        sigma_lux[mass][m.name] = lux.sigma_limit(mass=mass, sigma_guess=sigma_guess[m.name],
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[lux.name])
        sigma_pandax[mass][m.name] = lux.sigma_limit(mass=mass, sigma_guess=sigma_guess[m.name],
                                                  sigma_name=sigma_names[m.name],
                                                  fnfp_name=fnfp_names[m.name],
                                                  fnfp_val=fnfp_vals[m.name],
                                                  Nbackground=Nbg[pandax.name])
        sigma_cdmslite[mass][m.name] = cdmslite.sigma_limit(mass=mass, sigma_guess=sigma_guess[m.name],
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[cdmslite.name])
        sigma_supercdms[mass][m.name] = supercdms.sigma_limit(mass=mass, sigma_guess=sigma_guess[m.name],
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[supercdms.name])
        sigma_kimsiodLO[mass][m.name] = kimsiodLO.sigma_limit(mass=mass, sigma_guess=sigma_guess[m.name],
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[kimsiodLO.name])
        sigma_kimsiodHI[mass][m.name] = kimsiodHI.sigma_limit(mass=mass, sigma_guess=sigma_guess[m.name],
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[kimsiodHI.name])
        sigma_picoA[mass][m.name] = picoA.sigma_limit(mass=mass, sigma_guess=sigma_guess[m.name],
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[picoA.name])
        sigma_picoB[mass][m.name] = picoB.sigma_limit(mass=mass, sigma_guess=sigma_guess[m.name],
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[picoB.name])
        sigma_picoC[mass][m.name] = picoC.sigma_limit(mass=mass, sigma_guess=sigma_guess[m.name],
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[picoC.name])
        sigma_picoD[mass][m.name] = picoD.sigma_limit(mass=mass, sigma_guess=sigma_guess[m.name],
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[picoD.name])
        sigma_lims[m.name][i] = min(x for x in [sigma_supercdms[mass][m.name],
                                     sigma_cdmslite[mass][m.name],
                                     sigma_lux[mass][m.name],
                                     sigma_pandax[mass][m.name],
                                     sigma_kimsiodLO[mass][m.name],
                                     sigma_kimsiodHI[mass][m.name],
                                     sigma_picoA[mass][m.name],
                                     sigma_picoB[mass][m.name],
                                     sigma_picoC[mass][m.name],
                                     sigma_picoD[mass][m.name]] if x > 0.)

f_sigma_lims = {}
any_sigma_lims = {}
for m in MODELS_UV1 + MODELS_UV2:
    any_sigma_lims[m.name] = interpolate(np.log10(MASSES), np.log10(sigma_lims[m.name]),s=0)
def f_sigma_lims(modelname, mass):
    return 10**any_sigma_lims[modelname](np.log10(mass))



################################################################################################

def line_plots_1vsall(nsim=100, startsim=1, masses=[50.],
                      experiment_names=['Xe'],#['I','Ilo'], ['Xe','Xelo','Xehi','Xewide']
                      simmodels=[SI_Higgs], models=[SI_Higgs, anapole], time_info='True', GF=False,
                      filelabel='No_Time', allverbose=True, verbose=True,
                      results_root='/Users/SamWitte/Desktop/dmdd/Storage/results_uv/',
                      saveplots=True, alpha=0.3, xoffset=0.1, fs=20, fs2=18, sigma_lim_file=None, colors_list=None,
                      timeonly=False):

    if colors_list == None:
        for i,experiment in enumerate(experiments):
            colors_list = Colors[experiment_names[i]]

    if time_info == 'Both':
        time_list = ['T', 'F']
    elif time_info == 'True':
        time_list = ['T']
    elif time_info == 'False':
        time_list =['F']
    else:
        print 'ERROR'

    if len(filelabel)>0:
        filelabel = '_' + filelabel
    if sigma_lim_file is None:
        sigma_pandax = {}
        sigma_limvals = {}
        sigma_lux = {}
        sigma_supercdms = {}
        sigma_cdmslite = {}
    else:
        sigma_limvals = pickle.load(open(sigma_lim_file,'rb'))
    experiment_labels = []#experiment_names
    experiments = []
    for i,exn in enumerate(experiment_names):
        experiments.append(ALL_EXPERIMENTS[exn])
        experiment_labels.append(Experiment_LaTeX[exn])

    if time_info == 'Both':
        indexchange = 0
        for x in range(0, len(experiment_labels)):
            experiment_labels = np.append(experiment_labels,experiment_labels[x] + '\n No Time')
            indexchange += 1

    elif time_info == 'False':
        for x in range(0, len(experiment_labels)):
            experiment_labels[x] = experiment_labels[x] + '\n No Time'
    if time_info == 'Both':
        holdarray = copy.copy(experiment_labels)
        for x in range(0,experiment_labels.size / 2):
            holdarray[2*x]=experiment_labels[x]
            holdarray[2*x+1] = experiment_labels[x+experiment_labels.size/2]

        experiment_labels=holdarray

    for mass in masses:
        print '{}GeV:'.format(mass)

        if sigma_lim_file is None:
            sigma_limvals[mass] = {}
            sigma_lux[mass] = {}
            sigma_pandax[mass] = {}
            sigma_cdmslite[mass] = {}
            sigma_supercdms[mass] = {}
            for m in models:
                sigma_lux[mass][m.name] = lux.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                          fnfp_name=fnfp_names[m.name], fnfp_val=fnfp_vals[m.name],
                                                          Nbackground=Nbg[lux.name])
                sigma_pandax[mass][m.name] = pandax.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                          fnfp_name=fnfp_names[m.name], fnfp_val=fnfp_vals[m.name],
                                                          Nbackground=Nbg[pandax.name])
                sigma_cdmslite[mass][m.name] = cdmslite.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                                    fnfp_name=fnfp_names[m.name], fnfp_val=fnfp_vals[m.name],
                                                                    Nbackground=Nbg[cdmslite.name])
                sigma_supercdms[mass][m.name] = supercdms.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                                      fnfp_name=fnfp_names[m.name], fnfp_val=fnfp_vals[m.name],
                                                                      Nbackground=Nbg[supercdms.name])
                sigma_limvals[mass][m.name]=min(sigma_supercdms[mass][m.name],
                                                sigma_cdmslite[mass][m.name],
                                                sigma_lux[mass][m.name],
                                                sigma_pandax[mass][m.name])
                #   sigma_limvals[mass][m.name] = f_sigma_lims(m.name,mass)




        for m in simmodels:
            print ''
            print ''
            print m.name
            plt.figure()
            ax=plt.gca()
            for tval in time_list:
                if time_info == 'Both':
                    ticks = np.arange(2 * len(experiments)) + 0.5
                else:
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

                            mnrun = MultinestRun('sim{}'.format(n), experimentlist, m, pardic,
                                                    fitm, prior_ranges=prior_ranges,
                                                    force_sim=False,n_live_points=n_live_points,silent=True,
                                                    time_info=tval, GF=GF, TIMEONLY=timeonly)
                            if allverbose:
                                print ''
                                print mnrun.foldername
                            if verbose and not allverbose and fitm.name == m.name:
                                print ''
                                print '({})'.format(mnrun.foldername)

                            ln_ev[j] = mnrun.get_evidence()
                            renorm = np.amax(np.array([renorm,ln_ev[j]]))

                        for j,fitm in enumerate(models):
                            print j, ln_ev[j]
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

                    if time_info != 'Both':
                        for y in ys:
                            plt.plot( [ i + xoffset, i + 1 - xoffset ], [ y, y ], lw=1,
                                      alpha=alpha, color=colors_list[i])
                        plt.plot( [ i + 0.5, i + 0.5 ], [ ylo, yhi ], color='k', lw=2)
                        plt.plot( i + 0.5, ymedian, 'x', ms = 10, color='k', mew=2 )
                        success = float(np.sum( ys > 0.9 )) / len( ys ) * 100.

                        if success > 0:
                            plt.annotate( '{:.0f}\\%%'.format(success), xy = ( i + 0.5, 1.05 ),
                                          fontsize=fs2, va='center', ha='center', color='k' )
                    else:
                        if tval == 'F':
                            ii = 2*i + 1
                        else:
                            ii=2*i

                        for y in ys:
                            plt.plot( [ ii + xoffset, ii + 1 - xoffset ], [ y, y ], lw=1.25,
                                      alpha=alpha, color=colors_list[ii])
                        plt.plot( [ ii + 0.5, ii + 0.5 ], [ ylo, yhi ], color='k', lw=2)
                        plt.plot( ii + 0.5, ymedian, 'x', ms = 10, color='k', mew=2 )
                        success = float(np.sum( ys > 0.9 )) / len( ys ) * 100.

                        if success > 0:
                            plt.annotate( '{:.0f}\\%%'.format(success), xy = ( ii + 0.5, 1.05 ),
                                          fontsize=fs2, va='center', ha='center', color='k' )


            ax.set_xticks( ticks )
            ax.set_xticklabels( experiment_labels, size='small' )
            ax.set_title('True model: {} (mass: {:.0f} GeV)'.format(MODELNAME_TEX[m.name], mass), fontsize=fs)
            pl.ylim( ( -0.04, 1.19 ) )
            pl.axhline(0.9, ls='--', color='k', lw=2)

            pl.annotate('Simulations preferring true model (with \\textgreater 0.9 probability):',
                        xy=( 0.02, 0.96 ), xycoords='axes fraction', fontsize=fs2, va='top', color='k')
            pl.ylabel('Probability of true model', fontsize=fs)

            figname = results_root + 'lineplot_{:.0f}GeV_{}_{}sims{}.pdf'.format(mass, m.name, nsim, filelabel)
            if saveplots:
                pl.savefig(figname)



#########################################
def OneDhistogram(nsim=50, startsim=1, masses=[50.],
                  experiment_names=['Xe'],#['I','Ilo'], ['Xe','Xelo','Xehi','Xewide']
                  simmodels=[SI_Higgs], models=[SI_Higgs, anapole], time_info='Both', GF=True,
                  hspace = (1.06 * 50. ** (-1. / 5.)), filelabel='', allverbose=True, verbose=True,
                  results_root=os.environ['DMDD_AM_MAIN_PATH']+'/results_uv/', timeonly=False,
                  saveplots=True, alpha=0.3, xoffset=0.1, fs=18, fs2=18, sigma_lim_file=None,
                  colors_list=['Blue','Red','Green','Magenta','Orange']):

    xlinspace = np.linspace(0, 1, 200)
    leg_top = 0.
    leg_down = 0.
    success_list = []
    label_list = []
    ymax_list = np.zeros(4)
    leg_down = np.zeros(4)
    leg_top = np.zeros(4)
    lab = []

    if colors_list == None:
        for i,experiment in enumerate(experiments):
            colors_list = Colors[experiment_names[i]]

    totlines = len(experiment_names)
    if time_info == 'Both':
        time_list = ['F', 'T']
        totlines *= 2
    elif time_info == 'True':
        time_list = ['T']
    elif time_info == 'False':
        time_list =['F']
    else:
        print 'ERROR'


    if len(filelabel) > 0:
        filelabel = '_' + filelabel
    if sigma_lim_file is None:
        sigma_limvals = {}
        sigma_lux = {}
        sigma_supercdms = {}
        sigma_cdmslite = {}
        sigma_pandax = {}
    else:
        sigma_limvals = pickle.load(open(sigma_lim_file,'rb'))
    experiment_labels = []#experiment_names
    experiments = []
    for i,exn in enumerate(experiment_names):
        experiments.append(ALL_EXPERIMENTS[exn])
        experiment_labels.append(Experiment_LaTeX[exn])

    if time_info == 'Both':
        indexchange = 0
        for x in range(0, len(experiment_labels)):
            experiment_labels = np.append(experiment_labels, experiment_labels[x] + '\n No Time')
            indexchange += 1

    elif time_info == 'False':
        for x in range(0, len(experiment_labels)):
            experiment_labels[x] = experiment_labels[x] + '\n No Time'
    if time_info == 'Both':
        holdarray = copy.copy(experiment_labels)
        for x in range(0, len(experiment_labels) / 2):
            holdarray[2 * x] = experiment_labels[x]
            holdarray[2 * x + 1] = experiment_labels[x + experiment_labels.size / 2]

        experiment_labels = holdarray

    for mass in masses:
        print '{}GeV:'.format(mass)

        if sigma_lim_file is None:
            sigma_limvals[mass] = {}
            sigma_lux[mass] = {}
            sigma_cdmslite[mass] = {}
            sigma_supercdms[mass] = {}
            sigma_pandax[mass] = {}
            for m in models:
                sigma_lux[mass][m.name] = lux.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                          fnfp_name=fnfp_names[m.name], fnfp_val=fnfp_vals[m.name],
                                                          Nbackground=Nbg[lux.name])
                sigma_pandax[mass][m.name] = pandax.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                          fnfp_name=fnfp_names[m.name], fnfp_val=fnfp_vals[m.name],
                                                             Nbackground=Nbg[pandax.name])
                sigma_cdmslite[mass][m.name] = cdmslite.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                                    fnfp_name=fnfp_names[m.name], fnfp_val=fnfp_vals[m.name],
                                                                    Nbackground=Nbg[cdmslite.name])
                sigma_supercdms[mass][m.name] = supercdms.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                                      fnfp_name=fnfp_names[m.name], fnfp_val=fnfp_vals[m.name],
                                                                      Nbackground=Nbg[supercdms.name])
                sigma_limvals[mass][m.name]=min(sigma_supercdms[mass][m.name],
                                                sigma_cdmslite[mass][m.name],
                                                sigma_lux[mass][m.name],
                                                sigma_pandax[mass][m.name])
                #   sigma_limvals[mass][m.name] = f_sigma_lims(m.name,mass)


        for m in simmodels:
            print ''
            print ''
            print m.name
            fig = plt.figure()
            f, axarr = plt.subplots(2, 2, sharex='col')
            for tval in time_list:
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

                            mnrun = MultinestRun('sim{}'.format(n), experimentlist, m, pardic,
                                                    fitm, prior_ranges=prior_ranges,
                                                    force_sim=False,n_live_points=n_live_points,silent=True,
                                                    time_info=tval, GF=GF, TIMEONLY=timeonly)
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

                    if tval == 'F':
                        ii = 2*i + 1
                    else:
                        ii=2*i

                    probdistr = np.zeros_like(xlinspace)
                    std_dev = np.std(ys)
                    hspace = 1.06 * std_dev * nsim ** (-1. / 5.)
                    for x in range(0, xlinspace.size):
                        probdistr[x] = (1. / (nsim * hspace)) * norm.pdf((xlinspace[x] - ys) / hspace).sum()

                    success = float(np.sum(ys > 0.9)) / len(ys) * 100.
                    probdistr = probdistr
                    #success_list.append(success)
                    #label = experiment_labels[ii]
                    #label_list.append(experiment_labels[ii])
                    #c_list.append(colors_list[ii])
                    if i == 0:
                        ax_x = 0
                        ax_y = 0
                    elif i == 1:
                        ax_x = 0
                        ax_y = 1
                    elif i == 2:
                        ax_x = 1
                        ax_y = 0
                    else:
                        ax_x = 1
                        ax_y = 1
                    ii = i
                    if tval == 'F':
                        axarr[ax_x, ax_y].set_xlim([10., 100.])
                        ticks = np.power(10., 2 * np.array([.5, .7, .8, .9, .95]))
                        axarr[ax_x, ax_y].axes.get_xaxis().set_ticks(ticks)
                        axarr[ax_x, ax_y].axes.get_xaxis().set_ticklabels(['{:.0f}'.format(50),
                                                                           '{:.0f}'.format(70),
                                                                           '{:.0f}'.format(80),
                                                                           '{:.0f}'.format(90),
                                                                           '{:.0f}'.format(95)])
                        ls = '--'
                        ymax_list[i] = np.max(probdistr) + 0.1
                        axarr[ax_x, ax_y].plot(10 ** (2. * xlinspace), probdistr, ls, linewidth=2,
                                               color=colors_list[ii], dashes=(10, 10))
                        lab.append('(Dashed) No Time: ' + r'[Success: {:.0f}$\%$]'.format(success))
                        axarr[ax_x, ax_y].axes.get_yaxis().set_ticks([])
                    else:
                        ymax_list[i] = np.max([ymax_list[i], np.max(probdistr) + 0.1])

                        axarr[ax_x, ax_y].plot(10 ** (2. * xlinspace), probdistr, linewidth=2, color=colors_list[ii])
                        lab[i] += ('\n (Solid) Time: ' + r'[Success: {:.0f}$\%$]'.format(success))
                        if i > 1:
                            leg_top[i] = 0.7 * ymax_list[i]
                        else:
                            leg_top[i] = 0.1 * ymax_list[i]
                        axarr[ax_x, ax_y].text(12, leg_top[i], lab[i], color=colors_list[i], fontsize=10)
                        axarr[ax_x, ax_y].text(12, 0.85 * ymax_list[i], experiment_labels[2 * ii], color='k', fontsize=16)
                        axarr[ax_x, ax_y].set_ylim([0., ymax_list[i]])

            #axarr[1, 1].text(10, 0.5 * ymax, 'Solid: Time \n Dashed: No Time', color='k', fontsize=10)
            pl.suptitle('True model: {} (mass: {:.0f} GeV)'.format(MODELNAME_TEX[m.name], mass), fontsize=fs)
            f.text(0.5, .05, r'Probability of True Model   [$\%$]', ha='center', va='center', fontsize=fs)
            f.tight_layout(rect=(0, .05, 1, .95))
            # fig.text(0.06, 0.5, 'common ylabel', ha='center', va='center', rotation='vertical')
            plt.subplots_adjust(wspace=0.0, hspace=0.)
            figname = results_root + 'PDF_{:.0f}GeV_{}_{}sims{}.pdf'.format(mass, m.name, nsim, filelabel)
            if saveplots:
                pl.savefig(figname)




#########################################
def OneDhistogram_timeDiff(nsim=50, startsim=1, masses=[50.],
                  experiment_names=['Xe'],#['I','Ilo'], ['Xe','Xelo','Xehi','Xewide']
                  simmodels=[SI_Higgs], models=[SI_Higgs, anapole], GF=True,
                  hspace = (1.06 * 50. ** (-1. / 5.)), filelabel='', allverbose=True, verbose=True,
                  results_root=os.environ['DMDD_AM_MAIN_PATH']+'/results_uv/', timeonly=False,
                  saveplots=True, alpha=0.3, xoffset=0.1, fs=16, fs2=18, sigma_lim_file=None,
                  colors_list=['Aqua','Red','Black','Green','Magenta','Orange']):

    xlinspace = np.linspace(-.5,.5,300)
    leg_top = 0.
    leg_down = 0.
    avg_list = []
    med_list = []
    label_list = []
    c_list = []

    if colors_list == None:
        for i,experiment in enumerate(experiments):
            colors_list = Colors[experiment_names[i]]


    if len(filelabel)>0:
        filelabel = '_' + filelabel
    if sigma_lim_file is None:
        sigma_limvals = {}
        sigma_lux = {}
        sigma_supercdms = {}
        sigma_cdmslite = {}
        sigma_pandax = {}
    else:
        sigma_limvals = pickle.load(open(sigma_lim_file,'rb'))
    experiment_labels = []#experiment_names
    experiments = []
    for i,exn in enumerate(experiment_names):
        experiments.append(ALL_EXPERIMENTS[exn])
        experiment_labels.append(Experiment_LaTeX[exn])


    for mass in masses:
        print '{}GeV:'.format(mass)

        if sigma_lim_file is None:
            sigma_limvals[mass] = {}
            sigma_lux[mass] = {}
            sigma_cdmslite[mass] = {}
            sigma_supercdms[mass] = {}
            sigma_pandax[mass] = {}
            for m in models:
                sigma_lux[mass][m.name] = lux.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                          fnfp_name=fnfp_names[m.name], fnfp_val=fnfp_vals[m.name],
                                                          Nbackground=Nbg[lux.name])
                sigma_pandax[mass][m.name] = pandax.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                                fnfp_name=fnfp_names[m.name], fnfp_val=fnfp_vals[m.name],
                                                                Nbackground=Nbg[pandax.name])
                sigma_cdmslite[mass][m.name] = cdmslite.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                                    fnfp_name=fnfp_names[m.name], fnfp_val=fnfp_vals[m.name],
                                                                    Nbackground=Nbg[cdmslite.name])
                sigma_supercdms[mass][m.name] = supercdms.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                                      fnfp_name=fnfp_names[m.name], fnfp_val=fnfp_vals[m.name],
                                                                      Nbackground=Nbg[supercdms.name])
                sigma_limvals[mass][m.name]=min(sigma_supercdms[mass][m.name],
                                                sigma_cdmslite[mass][m.name],
                                                sigma_lux[mass][m.name],
                                                sigma_pandax[mass][m.name])
                #   sigma_limvals[mass][m.name] = f_sigma_lims(m.name,mass)




        for m in simmodels:
            print ''
            print ''
            print m.name
            plt.figure()
            ax=plt.gca()

            for i,experiment in enumerate(experiments):
                try:
                    if len(experiment) > 1:
                        experimentlist = experiment
                except:
                    experimentlist = [experiment]
                ys = np.zeros(nsim)
                ys_nt = np.zeros(nsim)

                for ni,n in enumerate(np.arange(startsim, nsim + startsim)):
                    renorm = -np.inf
                    ev_sum = 0.
                    ln_ev = np.zeros(len(models))
                    prob = np.zeros(len(models))

                    renorm_nt = -np.inf
                    ev_sum_nt = 0.
                    ln_ev_nt = np.zeros(len(models))
                    prob_nt = np.zeros(len(models))

                    for j,fitm in enumerate(models):
                        pardic = {'mass': mass, sigma_names[m.name]: sigma_limvals[mass][m.name]}

                        mnrun = MultinestRun('sim{}'.format(n), experimentlist, m, pardic,
                                                fitm, prior_ranges=prior_ranges,
                                                force_sim=False,n_live_points=n_live_points,silent=True,
                                                time_info='T', GF=GF, TIMEONLY=timeonly)

                        mnrun_nt = MultinestRun('sim{}'.format(n), experimentlist, m, pardic,
                                                fitm, prior_ranges=prior_ranges,
                                                force_sim=False,n_live_points=n_live_points,silent=True,
                                                time_info='F', GF=GF, TIMEONLY=timeonly)

                        if allverbose:
                            print ''
                            print mnrun.foldername
                        if verbose and not allverbose and fitm.name == m.name:
                            print ''
                            print '({})'.format(mnrun.foldername)

                        ln_ev[j] = mnrun.get_evidence()
                        renorm = np.amax(np.array([renorm,ln_ev[j]]))

                        ln_ev_nt[j] = mnrun_nt.get_evidence()
                        renorm_nt = np.amax(np.array([renorm_nt,ln_ev_nt[j]]))

                    for j,fitm in enumerate(models):
                        ev_sum += np.exp(-renorm + ln_ev[j])
                        ev_sum_nt += np.exp(-renorm_nt + ln_ev_nt[j])

                    for j,fitm in enumerate(models):
                        prob[j] = np.exp(-renorm + ln_ev[j]) / ev_sum
                        prob_nt[j] = np.exp(-renorm_nt + ln_ev_nt[j]) / ev_sum_nt
                        if fitm.name != m.name:
                            if verbose or allverbose:
                                print 'Time ', '{} {}'.format(fitm.name, prob[j])
                                print 'No Time ', '{} {}'.format(fitm.name, prob_nt[j])
                        else:
                            ys[ni] = prob[j]
                            ys_nt[ni] = prob_nt[j]
                            if verbose or allverbose:
                                print 'Time ','{} {} <--sim'.format(fitm.name, prob[j])
                                print 'No Time ','{} {} <--sim'.format(fitm.name, prob_nt[j])
                            if verbose and prob[j]<0.1:
                                print 'FLAG! truth unlikely?   <-----------|'
                                print ''


                    if verbose:
                        print ''

                    if ev_sum == np.inf:
                        print 'yikes. evidence is infinity.'
                ys = ys - ys_nt
                #for y in ys:
                #    plt.plot([y,y], [ 0., 100. ], lw=1,
                #             alpha=alpha, color=colors_list[i])

                probdistr = np.zeros(xlinspace.size)
                std_dev = np.std(ys)
                hspace = 1.06 * std_dev * nsim ** (-1. / 5.)

                for x in range(0, xlinspace.size):
                    probdistr[x] = (1. / (nsim * hspace)) * norm.pdf((xlinspace[x] - ys) / hspace).sum()
                avg = np.mean(ys) * 100.
                med = np.median(ys) * 100.
                avg_list.append(avg)
                med_list.append(med)
                label = experiment_labels[i]
                label_list.append(experiment_labels[i])
                c_list.append(colors_list[i])

                label = experiment_labels[i]
                #plt.plot(xlinspace, probdistr, linewidth=1, color=colors_list[i], label=label)

                bins = np.linspace(-.1, .25, 8)
                plt.hist(ys, bins, ec=colors_list[i], histtype='step')

                maxval = np.max(probdistr) + 1
                if i == 0:
                    maxylim = np.max(probdistr) + 1
                else:
                    if maxval > maxylim:
                        maxylim = maxval

            leg_down = maxylim / 12.
            leg_top = maxylim - leg_down
            for i in range(len(label_list)):
                plt.text(0.3, leg_top, label_list[i] +
                         r'  [$\bar{{\Delta}}$ = {:.1f} \%, $\tilde{{\Delta}}$ = {:.1f} \%]'.format(avg_list[i],
                                                                                                med_list[i]),
                         color=c_list[i], fontsize=10)
                leg_top -= leg_down
            ax.set_title('True model: {} (mass: {:.0f} GeV)'.format(MODELNAME_TEX[m.name], mass), fontsize=fs)

            pl.xlim([-.05, .25])
            #pl.ylim([0., maxylim])
            ax.axes.get_yaxis().set_ticks([])

            pl.ylabel('Density', fontsize=fs)
            pl.xlabel('$\Delta$ Probability', fontsize=fs)

            figname = results_root + 'PDF_{:.0f}GeV_{}_{}sims{}.pdf'.format(mass, m.name, nsim, filelabel)
            if saveplots:
                pl.savefig(figname)




###################################

def compare_UV_models(masses=[50.], experiment_name='Xe',
                   models=MODELS_UV1, plottag='models1',
                   sigma_limvals = None, separate_experiments=False,
                   results_root='/Users/SamWitte/Desktop/dmdd/Storage/results_uv/',
                   xmin=1,xmax=100, ymax= None, fontsize = 20,
                   show_legend=True, save_plot=True,
                   v_lag=220, v_rms=220, v_esc=544, rho_x=0.3):

    linest = [ '-', '--','-.', ':','-', '--','-.', ':']
    linwid = [4,4,4,4,3,3,3,3]
    if type(ALL_EXPERIMENTS[experiment_name]) == list:
        experiments = ALL_EXPERIMENTS[experiment_name]
    else:
        experiments = [ALL_EXPERIMENTS[experiment_name]]
    qs = np.linspace(xmin,xmax,1000)
    plt.figure()

    for mi,mass in enumerate(masses):
        if sigma_limvals is None:
            sigma_limvals = {}

        for j,m in enumerate(models):
            if m.name not in sigma_limvals.keys():
                sigma_limvals[m.name] = pandax.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                fnfp_name=fnfp_names[m.name],
                                                fnfp_val=fnfp_vals[m.name],
                                                           Nbackground=Nbg[pandax.name])

            kwargs = {
                'mass': mass,
                sigma_names[m.name]: sigma_limvals[m.name],
                'v_lag': v_lag,
                'v_rms': v_rms,
                'v_esc': v_esc,
                'rho_x': rho_x,
                }
            if len(m.fixed_params) > 0:
                for k,v in m.fixed_params.iteritems():
                    kwargs[k] = v

            Nexp_tot = np.zeros(len(qs))

            for i,experiment in enumerate(experiments):
                allkwargs = dict(kwargs)
                allkwargs['element'] = experiment.element
                drdq = rate_UV.dRdQ( qs, **allkwargs )
                Nexp = drdq * experiment.exposure * YEAR_IN_S * experiment.efficiency( qs )
                Nexp_tot += Nexp
                if separate_experiments:
                    plt.plot(qs, Nexp, linest[i], color=MODELNAME_COLOR[m.name], lw=1)
            if mi==0:
                label = MODELNAME_TEX[m.name]
            else:
                label=None
            plt.plot(qs, Nexp, linestyle=linest[j], color=MODELNAME_COLOR[m.name], label=label, lw=linwid[j])

    ###

    plt.xlim(xmax=xmax)
    plt.xlim(xmin=xmin)
    plt.ylim(ymax=ymax)
    plt.ylim(ymin=0)

    xlabel = 'Nuclear recoil energy [keV]'
    ylabel = 'Number of events'
    ax = plt.gca()
    fig = plt.gcf()
    xlabel = ax.set_xlabel(xlabel,fontsize=fontsize)
    ylabel = ax.set_ylabel(ylabel,fontsize=fontsize)
    plt.title('Recoil spectra on {} (WIMP mass: {:.0f}GeV)'.format(experiment_name,mass), fontsize=fontsize)
    if show_legend:
        legend = ax.legend(prop={'size':20}, loc=[1,0], fontsize=fontsize)
    if save_plot:
        filename = results_root + 'all_{}_{:.0f}GeV'.format(plottag, mass)
        for experiment in experiments:
            filename += '_{}'.format(experiment.name)
        filename += '.pdf'
        if show_legend:
            pl.savefig(filename, bbox_extra_artists=[xlabel, ylabel, legend], bbox_inches='tight')
        else:
            pl.savefig(filename, bbox_extra_artists=[xlabel, ylabel], bbox_inches='tight')
        #plt.savefig(filename)





################################################################################################

def line_plots_classes(nsim=50, startsim=1, masses=[50.], experiment_names=['Xe'],
                      allmodels=MODELS_UV1, classes=MODEL_CLASSES, classes_names=MODEL_CLASSES_NAMES,
                      filelabel='', verbose=True,
                      results_root='/Users/SamWitte/Desktop/dmdd/Storage/results_uv/',
                      saveplots=True, alpha=0.3, xoffset=0.1, fs=20, fs2=18):

    if len(filelabel)>0:
        filelabel = '_' + filelabel
    sigma_limvals = {}
    sigma_lux = {}
    sigma_supercdms = {}
    sigma_cdmslite = {}
    experiment_labels = experiment_names
    experiments = []
    for exn in experiment_names:
        experiments.append(ALL_EXPERIMENTS[exn])

    for mass in masses:
        print '{}GeV:'.format(mass)
        sigma_limvals[mass] = {}
        sigma_lux[mass] = {}
        sigma_cdmslite[mass] = {}
        sigma_supercdms[mass] = {}
        for m in allmodels:
            sigma_lux[mass][m.name] = lux.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                      fnfp_name=fnfp_names[m.name], fnfp_val=fnfp_vals[m.name],
                                                      Nbackground=Nbg[lux.name])
            sigma_pandax[mass][m.name] = lux.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                        fnfp_name=fnfp_names[m.name], fnfp_val=fnfp_vals[m.name],
                                                         Nbackground=Nbg[pandax.name])
            sigma_cdmslite[mass][m.name] = cdmslite.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                                fnfp_name=fnfp_names[m.name], fnfp_val=fnfp_vals[m.name],
                                                                Nbackground=Nbg[cdmslite.name])
            sigma_supercdms[mass][m.name] = supercdms.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                                  fnfp_name=fnfp_names[m.name], fnfp_val=fnfp_vals[m.name],
                                                                  Nbackground=Nbg[supercdms.name])
            sigma_limvals[mass][m.name]=min(sigma_supercdms[mass][m.name],
                                            sigma_cdmslite[mass][m.name],
                                            sigma_lux[mass][m.name],
                                            sigma_pandax[mass][m.name])




        #sigma_limvals[mass] = {}
        #for m in allmodels:
        #    sigma_limvals[mass][m.name] = lux.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
        #                                                    fnfp_name=fnfp_names[m.name],
        #                                                    fnfp_val=fnfp_vals[m.name])

        for m in allmodels:
            print 'truth: ', m.name
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
                    for mi,mclass in enumerate(classes_names):
                        if m.name in mclass:
                            trueclass_names = mclass
                            trueclass = classes[mi]
                            trueclass_factor = 1.#/ len(trueclass_names) #float(len(allmodels)) / (len(trueclass_names) * len(classes))
                            #print 'true class: ',trueclass_factor, trueclass_names
                    renorm = -np.inf
                    ev_sum = 0.
                    ln_ev = np.zeros(len(allmodels))
                    prob = np.zeros(len(allmodels))

                    for j,fitm in enumerate(allmodels):
                        pardic = {'mass': mass, sigma_names[m.name]: sigma_limvals[mass][m.name]}

                        mnrun = MultinestRun('sim{}'.format(n), experimentlist, m, pardic,
                                                fitm, prior_ranges=prior_ranges,
                                                force_sim=False,n_live_points=n_live_points,silent=True)
                        ln_ev[j] = mnrun.get_evidence()
                        renorm = np.amax(np.array([renorm,ln_ev[j]]))

                    for j,fitm in enumerate(allmodels):
                        ev_sum += np.exp(-renorm + ln_ev[j])

                    for j,fitm in enumerate(allmodels):
                        prob[j] = np.exp(-renorm + ln_ev[j]) / ev_sum
                        if fitm.name in trueclass_names:
                            ys[ni] += prob[j] * trueclass_factor
                            if verbose:
                                print 'adding {} {}'.format(fitm.name, prob[j])
                    if verbose:
                        print 'sim{}: trueclass prob={}'.format(n, ys[ni])

                    if ev_sum == np.inf:
                        print 'yikes. evidence is infinity.'


                ymedian = np.median(ys)
                ylo, yhi = np.percentile(ys, 25), np.percentile(ys, 75)
                for y in ys:
                    plt.plot( [ i + xoffset, i + 1 - xoffset ], [ y, y ], lw=1,
                              alpha=alpha, color=Colors[experiment_labels[i]])
                plt.plot( [ i + 0.5, i + 0.5 ], [ ylo, yhi ], color='k', lw=2)
                plt.plot( i + 0.5, ymedian, 'x', ms = 10, color='k', mew=2 )
                success = float(np.sum( ys > 0.9 )) / len( ys ) * 100.

                if success > 0:
                    plt.annotate( '{:.0f}\\%%'.format(success), xy = ( i + 0.5, 1.05 ),
                                  fontsize=fs2, va='center', ha='center', color='k' )

            ax.set_xticks( ticks )
            ax.set_xticklabels( experiment_labels )
            ax.set_title('True model: {} (mass: {:.0f} GeV)'.format(MODELNAME_TEX[m.name], mass), fontsize=fs)
            pl.ylim(( -0.04, 1.19 ))
            pl.axhline(0.9, ls='--', color='k', lw=2)

            pl.annotate('Simulations preferring true class (with \\textgreater 0.9 probability):',
                        xy=( 0.02, 0.96 ), xycoords='axes fraction', fontsize=fs2, va='top', color='k')
            pl.ylabel('Probability of true class', fontsize=fs)

            figname = results_root + 'lineplot_class_{:.0f}GeV_{}_{}sims{}.pdf'.format(mass, m.name, nsim,filelabel)
            if saveplots:
                pl.savefig(figname)


################################################################################################

def line_plots_all(nsim=1, startsim=1, masses=[50.], experiment_names=['Xe'],
                      allmodels=MODELS_UV1, classes=MODEL_CLASSES, classes_names=MODEL_CLASSES_NAMES,
                      filelabel='', verbose=True,
                      results_root='/Users/SamWitte/Desktop/dmdd/Storage/results_uv/',
                      saveplots=True, alpha=0.3, xoffset=0.1, fs=20, fs2=18):

    if len(filelabel)>0:
        filelabel = '_' + filelabel
    sigma_limvals = {}
    experiment_labels = experiment_names
    experiments = []
    for exn in experiment_names:
        experiments.append(ALL_EXPERIMENTS[exn])

    for mass in masses:
        sigma_limvals[mass] = {}
        for m in allmodels:
            sigma_limvals[mass][m.name] = pandax.sigma_limit(mass=mass, sigma_name=sigma_names[m.name],
                                                            fnfp_name=fnfp_names[m.name],
                                                            fnfp_val=fnfp_vals[m.name],
                                                             Nbackground=Nbg[pandax.name])

        for m in allmodels:
            print 'truth: ', m.name
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
                ysclass = np.zeros(nsim)

                for ni,n in enumerate(np.arange(startsim, nsim + startsim)):
                    for mi,mclass in enumerate(classes_names):
                        if m.name in mclass:
                            trueclass_names = mclass
                            trueclass = classes[mi]
                            trueclass_factor = 1.
                    renorm = -np.inf
                    ev_sum = 0.
                    ln_ev = np.zeros(len(allmodels))
                    prob = np.zeros(len(allmodels))

                    for j,fitm in enumerate(allmodels):
                        pardic = {'mass': mass, sigma_names[m.name]: sigma_limvals[mass][m.name]}

                        mnrun = MultinestRun('sim{}'.format(n), experimentlist, m, pardic,
                                                fitm, prior_ranges=prior_ranges,
                                                force_sim=False,n_live_points=n_live_points,silent=True)
                        ln_ev[j] = mnrun.get_evidence()
                        renorm = np.amax(np.array([renorm,ln_ev[j]]))

                    for j,fitm in enumerate(allmodels):
                        ev_sum += np.exp(-renorm + ln_ev[j])

                    for j,fitm in enumerate(allmodels):
                        prob[j] = np.exp(-renorm + ln_ev[j]) / ev_sum
                        if fitm.name == m.name:
                             ys[ni] = prob[j]
                        if fitm.name in trueclass_names:
                            ysclass[ni] += prob[j] * trueclass_factor
                            if verbose:
                                print 'adding {} {}'.format(fitm.name, prob[j])
                    if verbose:
                        print 'sim{}: trueclass prob={}, truemodel prob={} <--sim'.format(n, ysclass[ni], ys[ni])

                    if ev_sum == np.inf:
                        print 'yikes. evidence is infinity.'


                ymedian = np.median(ys)
                yclassmedian = np.median(ysclass)
                ylo, yhi = np.percentile(ys, 25), np.percentile(ys, 75)
                yclasslo, yclasshi = np.percentile(ysclass, 25), np.percentile(ysclass, 75)
                for y in ys:
                    plt.plot( [ i + xoffset, i + 1 - xoffset ], [ y, y ], lw=1,
                              alpha=alpha, color=Colors[experiment_labels[i]])
                plt.plot( [ i + 0.5, i + 0.5 ], [ ylo, yhi ], color='k', lw=2)
                plt.plot( i + 0.5, ymedian, 'x', ms = 10, color='k', mew=2 )

                plt.plot( [ i + 0.53, i + 0.53 ], [ yclasslo, yclasshi ], color='gray', lw=2,alpha=0.5)
                plt.plot( i + 0.53, yclassmedian, 'o', ms = 10, color='white', mew=2 )

                success = float(np.sum( ys > 0.9 )) / len( ys ) * 100.
                successclass = float(np.sum( ysclass > 0.9 )) / len( ysclass ) * 100.
                print successclass
                #if success > 0:
                plt.annotate( '{:.0f}\\%%'.format(success), xy = ( i + 0.5, 1.05 ),
                                  fontsize=fs2, va='center', ha='center', color='k' )
                plt.annotate( '({:.0f}\\%% )'.format(successclass), xy = ( i + 0.58, 1.05 ),
                                  fontsize=fs2, va='center', ha='center', color='gray' )

            ax.set_xticks( ticks )
            ax.set_xticklabels( experiment_labels )
            ax.set_title('True model: {} (mass: {:.0f} GeV)'.format(MODELNAME_TEX[m.name], mass), fontsize=fs)
            pl.ylim( ( -0.04, 1.19 ) )
            pl.axhline(0.9, ls='--', color='k', lw=2)

            pl.annotate('Simulations preferring true model (class):',
                        xy=( 0.02, 0.96 ), xycoords='axes fraction', fontsize=fs2, va='top', color='k')
            pl.ylabel('Probability of true model', fontsize=fs)

            figname = results_root + 'lineplot_classall_{:.0f}GeV_{}_{}sims{}.pdf'.format(mass, m.name, nsim,filelabel)
            if saveplots:
                pl.savefig(figname)




########################

folder1 = '/Users/SamWitte/Desktop/dmdd/Storage/chains_uv/sim1_Xe_fnfp_si_1.00_mass_50.00_sigma_si_67.91_fitfor_mass_sigma_sd_fixed_fnfp_sd_-1.10_logflat_nlive2000'
folder2 = '/Users/SamWitte/Desktop/dmdd/Storage/chains_uv/sim2_Xe_fnfp_si_1.00_mass_50.00_sigma_si_67.91_fitfor_mass_sigma_sd_fixed_fnfp_sd_-1.10_logflat_nlive2000'
folder3 = '/Users/SamWitte/Desktop/dmdd/Storage/chains_uv/sim3_Xe_fnfp_si_1.00_mass_50.00_sigma_si_67.91_fitfor_mass_sigma_sd_fixed_fnfp_sd_-1.10_logflat_nlive2000'

color_groups = [('Aqua', 'Blue'),
                ('Orange', 'Red'),
                ('Magenta', 'Pink')]

def plot_multiple_posteriors(folders=[folder1,folder3],
                             colors=None,
                             savefile=None, **kwargs):
    files = [f + '/1-post_equal_weights.dat' for f in folders]
    plt.figure()
    ax = plt.gca()
    for i,f in enumerate(files):
        plot_posterior(f, ax=ax, alpha=1,
                       contour_colors=color_groups[i],
                       zorder=10*i,
                       **kwargs)

    if savefile is not None:
        plt.savefig('test.png')

def plot_posterior(filename, **kwargs):
    x,y = np.loadtxt(filename, usecols=(0,1), unpack=True)
    return plot_2d_posterior(x,y, **kwargs)

def plot_2d_posterior(x, y,xlabel='', ylabel='',
                      input_x=None, input_y=None, input_color='red', fontsize=22,
                      contour_colors=('Aqua','Blue'), alpha=0.5,
                      xmin=None, xmax=None, ymin=None, ymax=None,
                      title='', plot_samples=False, samples_color='gray',
                      contour_lw=2, savefile=None, plot_contours=True,
                      show_legend=False, ax=None, return_ax=False,
                      zorder=0):
    """
    The main chains plotting routine, for visualizing 2d posteriors...
    """
    if ax is None:
        plt.figure()
        ax = pl.gca()

    n = 100

    points = np.array([x,y])
    posterior = kde(points)

    if xmin is None:
        xmin=x.min()

    if xmax is None:
        xmax=x.max()

    if ymin is None:
        ymin=y.min()

    if ymax is None:
        ymax=y.max()


    step_x = ( xmax - xmin ) / n
    step_y = ( ymax - ymin ) / n
    grid_pars = np.mgrid[0:n,0:n]
    par_x = grid_pars[0]*step_x + xmin
    par_y = grid_pars[1]*step_y + ymin
    grid_posterior = grid_pars[0]*0.

    for i in range(n):
        for j in range(n):
            grid_posterior[i][j] = posterior([par_x[i][j],par_y[i][j]])

    ax.set_title(title, fontsize=fontsize)
    ax.set_xlabel(xlabel,fontsize=fontsize)
    ax.set_ylabel(ylabel,fontsize=fontsize)
    if plot_samples:
        ax.plot(x,y,'o',ms=1, mfc=samplescolor, mec=samplescolor)
    ax.plot(input_x,input_y,'x',mew=3,ms=15,color=input_color,label='input',zorder=4)
    if plot_contours:
        percentage_integral = np.array([0.95,0.68,0.])
        contours = 0.* percentage_integral
        num_epsilon_steps = 1000.
        epsilon = grid_posterior.max()/num_epsilon_steps
        epsilon_marks = np.arange(num_epsilon_steps + 1)
        posterior_marks = grid_posterior.max() - epsilon_marks * epsilon
        posterior_norm = grid_posterior.sum()
        for j in np.arange(len(percentage_integral)):
            for i in epsilon_marks:
                posterior_integral = grid_posterior[np.where(grid_posterior>posterior_marks[i])].sum()/posterior_norm
                if posterior_integral > percentage_integral[j]:
                    break
            contours[j] = posterior_marks[i]
        contours[-1]=grid_posterior.max()
        ax.contour(par_x, par_y, grid_posterior, contours, linewidths=contour_lw,colors='k',zorder=3+zorder)
        ax.contourf(par_x, par_y, grid_posterior, contours,colors=contour_colors,alpha=alpha,zorder=2+zorder)

    ax.set_xlim(min(ax.get_xlim()[0],xmin), max(ax.get_xlim()[1],xmax))
    ax.set_xlim(min(ax.get_ylim()[0],ymin), max(ax.get_ylim()[1],ymax))
    if return_ax:
        return ax
    if show_legend:
        ax.legend(prop={'size':20},numpoints=1)
    if savefile is None:
        return par_x, par_y, grid_posterior, contours
    else:
        pl.savefig(savefile)


#############
############
##########
def make_Nexpected_latextable(masses=[15,50,500],experiments=[xe,ge,iod,flu],filename=None,
                              models=MODELS_UV1+MODELS_UV2,
                              results_root='/Users/SamWitte/Desktop/dmdd/Storage/results_uv/'):
    if filename is None:
        filename = results_root + 'latextable_Nexpected.txt'
    fout = open(filename, 'w')
    fout.write('\\begin{table*}[t]\n')
    fout.write('\\setlength{\\extrarowheight}{3pt}\n')
    fout.write('\\setlength{\\tabcolsep}{12pt}\n')
    fout.write('\\begin{center}\n')
    fout.write('\\begin{tabular}{c|')
    for experiment in experiments:
        fout.write('m{3cm}')
    fout.write('}\n')
    fout.write('Interaction [$\sigma_p$] /target')
    for experiment in experiments:
        fout.write('& {}'.format(experiment.name))
    fout.write('\\\\\n')
    fout.write('\\hline\\hline\\\\\n')

    for m in models:

        fout.write(MODELNAME_TEX[m.name])
        for experiment in experiments:
            entry = '& ('
            for mass in masses:
                Nexps = dmdd.Nexpected(experiment.element,
                                            experiment.Qmin,
                                            experiment.Qmax,
                                            experiment.exposure,
                                            experiment.efficiency,
                                            experiment.Start, experiment.End,
                                            sigma_names[m.name],
                                            f_sigma_lims(m.name,mass),
                                            mass=mass,
                                            fnfp_name=fnfp_names[m.name],
                                            fnfp_val=fnfp_vals[m.name])
                if mass == masses[0]:
                    entry += '{:.0f}'.format( Nexps )
                else:
                    entry += ', {:.0f}'.format( Nexps )
            entry += ')'
            fout.write(entry)
        fout.write('\\\\\n')

    fout.write('\\end{tabular}\n')
    fout.write('\\end{center}\n')
    fout.write('\\caption{}\n')
    fout.write('\\label{}\n')
    fout.write('\\end{table*}\n')

    fout.close()

##########################################
def get_chains(sim_name, experiments, sim_model,
            param_values, fit_model,
            prior='logflat',
            sim_root=SIM_PATH, chains_root=CHAINS_PATH,
            n_live_points=2000, basename='1-', time_tag='With_Time',
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

    foldername += '_{}_nlive{}_{}'.format(prior, n_live_points, time_tag)

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

def marginal_plot(samples, ax=None, orientation='vertical', **kwargs):

    n = 1000
    posterior = gaussian_kde(samples)
    minval=0.
    maxval=samples.max()
    step = (maxval-minval)/float(n)
    grid_par = np.arange(n)
    par = grid_par*step + minval
    grid_posterior = posterior(par)
    x = par
    y = grid_posterior

    if orientation=='horizontal':
        holder=x
        x=y*1.
        y=holder*1.

    if ax is None:
        ax = plt.gca()
    ax.plot(x, y, **kwargs)


def spaghetti(param, sim_range,mass=50,
            experiments=[xe,ge,flu], sim_model=SI_Higgs,
            fit_model=SI_Higgs,
            prior='logflat',xy=( 0.80, 0.90 ),labelfont=20,
            sim_root=SIM_PATH, chains_root=CHAINS_PATH,
            n_live_points=2000, basename='1-', time_tag='With_Time',
            show_xlabel=False, color='k',xmax=100,
            **kwargs):
    """

    param: name of parameter you want to plot (e.g., 'mass')
    sim_range: range of simulations to plot (e.g., [1,50])

    everything else is the same arguments with which you'd
    instantiate a MultinestRun

    additional keyword arguments passed to plt.plot
    """
    param_values = {}
    param_values['mass'] = mass
    param_values[sigma_names[sim_model.name]] = f_sigma_lims(sim_model.name,mass)
    plt.figure()
    ax = plt.gca()
    for i in range(*sim_range):
        sim_name = 'sim{}'.format(i)
        chain = get_chains(sim_name, experiments, sim_model,
                            param_values, fit_model,
                            prior=prior,sim_root=sim_root, chains_root=chains_root,
                            n_live_points=n_live_points, basename=basename, time_tag=time_tag,
                            params_to_get=param)

        marginal_plot(chain, ax=ax, color=color, **kwargs)
        ax.axvline(param_values[param], color='b')
        plt.xlim(xmax=xmax)

    plt.yticks([])
    #plt.ylabel('posterior probability density')
    pl.annotate(MODELNAME_TEX[sim_model.name], xy=xy,
                xycoords='axes fraction', fontsize=labelfont,
                va='top', color='k')
    if show_xlabel:
        plt.xlabel(PARAM_TEX[param], fontsize=22)

def all_spaghetti(models=MODELS_UV1,sim_range=[1,10],
                  tag='xegef',
                  param='mass',mass=50,
                  experiment=[xe,ge,flu],
                  saveplot=True,
                  prior='logflat',show_xlabel=False,color='k',
                  results_root='/Users/SamWitte/Desktop/dmdd/Storage/results_uv/',
                  xmax=None,
                  **kwargs):
    for m in models:
        spaghetti(param, sim_range,mass=50,
            experiments=experiment, sim_model=m,
            fit_model=m,
            prior=prior,
            sim_root=SIM_PATH, chains_root=CHAINS_PATH,
            n_live_points=2000, basename='1-',
            time_tag=time_tag, xmax=xmax,
            show_xlabel=show_xlabel, color=color,
            **kwargs)

        filename=results_root + 'spaghetti_{}_{}_{:.0f}GeV_{}.pdf'.format(param,m.name,mass,tag)
        plt.savefig(filename)


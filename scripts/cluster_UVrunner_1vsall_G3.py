#!/usr/bin/env python
#NOTE: this script expects that sims already exist. if some happen to not be there, it will create them; however, if all are absent, then there will be file-creation conflicts for sims the way this runner is set up.

import subprocess as sp
import os,sys,fnmatch
import argparse
import pickle 

import dmdd
import eff

import experiment_parameters as ep

parser = argparse.ArgumentParser()
parser.add_argument('--nofit', action='store_true')
parser.add_argument('--dosim', action='store_true')
parser.add_argument('--dumplims', action='store_true')
parser.add_argument('--limfile',default=None)
parser.add_argument('--tag',default='2G3+')
parser.add_argument('--masses',nargs='+',default=[50],type=float)
parser.add_argument('--ngroups', type=int, default=512)
parser.add_argument('--nsimgroups', type=int, default=512)
parser.add_argument('--nsim', type=int, default=1)
parser.add_argument('--startsim', type=int, default=1)
parser.add_argument('--nfit', type=int, default=1)
parser.add_argument('--startfit', type=int, default=1)
parser.add_argument('--prior',default='logflat')
parser.add_argument('--experiments',nargs='+',default=['XeG3','XeG3 I+ F+', 'I+ F+'])
parser.add_argument('--path',default='/home/verag/Projects/Repositories/dmdd_2014/scripts/')
#'/Users/verag/Research/Repositories/dmdd_2014/scripts/' #macbook


args = parser.parse_args()
DO_FIT = not args.nofit
DO_SIM = args.dosim
TAG = args.tag
MASSES = args.masses
NGROUPS = args.ngroups
NSIM = args.nsim
STARTSIM = args.startsim
NFIT = args.nfit
STARTFIT = args.startfit
PRIOR_TYPE = args.prior
NSIMGROUPS = args.nsimgroups 
SCRIPTS_PATH = args.path
EXPERIMENTS = args.experiments
SINGLE_EXPERIMENTS = []
for i,experiment in enumerate(EXPERIMENTS):
    labels = experiment.split()
    for lab in labels:
        if lab not in SINGLE_EXPERIMENTS:
            SINGLE_EXPERIMENTS.append(lab)
SINGLE_EXPERIMENTS.sort()
print SINGLE_EXPERIMENTS
DUMP_SIGMA_LIMS = args.dumplims
LIMFILE = args.limfile
            

SI_Higgs = dmdd.Haxton_Model('SI_Higgs', ['mass', 'sigma_si'], fixed_params={'fnfp_si': 1})
millicharge = dmdd.Haxton_Model('Millicharge', ['mass', 'sigma_si_massless'], fixed_params={'fnfp_si_massless': 0})
SD_flavoruniversal = dmdd.Haxton_Model('SD_fu', ['mass','sigma_sd'], fixed_params={'fnfp_sd': -1.1})
SD_Zmediated = dmdd.Haxton_Model('SD_Z-mediated', ['mass','sigma_sd'], fixed_params={'fnfp_sd': -13.2})
SD_Moira = dmdd.Haxton_Model('SI_Moira', ['mass','sigma_sd'], fixed_params={'fnfp_sd': 0.})
anapole = dmdd.Haxton_Model('Anapole', ['mass','sigma_anapole'])
magdip_heavy = dmdd.Haxton_Model('Mag.dip.heavy', ['mass','sigma_magdip'])
magdip_0 = dmdd.Haxton_Model('Mag.dip.light', ['mass','sigma_magdip_massless'])
elecdip_heavy = dmdd.Haxton_Model('Elec.dip.heavy', ['mass','sigma_elecdip'])
elecdip_0 = dmdd.Haxton_Model('Elec.dip.light', ['mass','sigma_elecdip_massless'])
f1 = dmdd.Haxton_Model('f1', ['mass','sigma_f1'], fixed_params={'fnfp_f1': 1.})
f2_Higgs = dmdd.Haxton_Model('f2_Higgs', ['mass','sigma_f2'], fixed_params={'fnfp_f2': -0.05})
f2_flavoruniversal = dmdd.Haxton_Model('f2_flavor-universal', ['mass','sigma_f2'], fixed_params={'fnfp_f2': 1.})
f3_Higgs = dmdd.Haxton_Model('f3_Higgs', ['mass','sigma_f3'], fixed_params={'fnfp_f3': -0.05})
f3_flavoruniversal = dmdd.Haxton_Model('f3_flavor-universal', ['mass','sigma_f3'], fixed_params={'fnfp_f3': 1.})
LS = dmdd.Haxton_Model('LS', ['mass','sigma_LS'], fixed_params={'fnfp_LS': 0.})



ALLMODELS = [SI_Higgs, millicharge, SD_flavoruniversal, SD_Zmediated, SD_Moira, anapole, magdip_heavy, magdip_0, elecdip_heavy, elecdip_0, f1, f2_Higgs, f2_flavoruniversal, f3_Higgs, f3_flavoruniversal]

MODELS1 = [SI_Higgs, millicharge, SD_flavoruniversal, anapole, magdip_heavy, magdip_0, elecdip_heavy, elecdip_0]
MODELS2 = [magdip_heavy, elecdip_heavy, f1, f2_Higgs, f2_flavoruniversal, f3_Higgs, f3_flavoruniversal, LS]

SIMMODELS = [f1, f2_Higgs, f2_flavoruniversal] #[SI_Higgs,magdip_heavy,magdip_0] 
FITMODELS = MODELS2


##get upper limits for a given mass, from G2 + current F and I:
experiment='Xe'
xe=dmdd.Experiment(experiment,ep.Target[experiment], ep.Qmin[experiment], ep.Qmax[experiment], ep.Exposure[experiment], eff.efficiency_Xe)
experiment = 'Ge'
ge=dmdd.Experiment(experiment,ep.Target[experiment], ep.Qmin[experiment], ep.Qmax[experiment], ep.Exposure[experiment], eff.efficiency_Xe)
kimsiodLO=dmdd.Experiment('KIMSiodLO','iodine',28.57,42.86, 32.8, eff.efficiency_KIMS)#Q=0.07
kimsiodHI=dmdd.Experiment('KIMSiodHI','iodine',42.86,300, 32.8, eff.efficiency_KIMS)#Q=0.07
picoA=dmdd.Experiment('PICO2LA','fluorine', 3.2, 400, 0.205, eff.efficiency_Xe)
picoB=dmdd.Experiment('PICO2LB','fluorine', 4.4, 400, 0.046, eff.efficiency_Xe)
picoC=dmdd.Experiment('PICO2LC','fluorine', 6.1, 400, 0.225, eff.efficiency_Xe)
picoD=dmdd.Experiment('PICO2LD','fluorine', 8.1, 400, 0.104, eff.efficiency_Xe)

current_experiments = [xe,ge,kimsiodLO,kimsiodHI,picoA,picoB,picoC,picoD]

Nbg = {
    xe.name: 1,
    ge.name: 1,
    kimsiodLO.name: 22,#500,
    kimsiodHI.name: 63,#4000,
    picoA.name: 9,
    picoB.name: 1,
    picoC.name: 3,
    picoD.name: 1,
}

mx_guess = {
    xe.name: 5.,
    ge.name: 1.,
    kimsiodLO.name: 1.,
    kimsiodHI.name: 1.,
    picoA.name: 1.,
    picoB.name: 1.,
    picoC.name: 1.,
    picoD.name: 1.,
}

sigma_guess = {}
for m in SIMMODELS + FITMODELS:
    sigma_guess[m.name] = 1e10
    #if m.name=='f2_Higgs' or m.name=='f3_Higgs':
    #    sigma_guess[m.name] = 1e6
    #else:
    #    sigma_guess[m.name] = 100.



sigma_names = {}
fnfp_names = {}
fnfp_vals = {}
for m in SIMMODELS + FITMODELS:
    sigma_names[m.name] = m.param_names[1]    
    if len(m.fixed_params)>0:
        fnfp_names[m.name] = m.fixed_params.keys()[0]
        fnfp_vals[m.name] = m.fixed_params.values()[0] 
    else:
        fnfp_names[m.name] = None
        fnfp_vals[m.name] = None
    
 
sigma_vals = {}
sigma_xe={}
sigma_ge={}
sigma_kimsiodLO={}
sigma_kimsiodHI={}
sigma_picoA={}
sigma_picoB={}
sigma_picoC={}
sigma_picoD={}

for i,mass in enumerate(MASSES):
    sigma_xe[mass]={}
    sigma_ge[mass]={}
    sigma_kimsiodLO[mass]={}
    sigma_kimsiodHI[mass]={}
    sigma_picoA[mass]={}
    sigma_picoB[mass]={}
    sigma_picoC[mass]={}
    sigma_picoD[mass]={}
    sigma_vals[mass] = {}
    for m in SIMMODELS:
        sigma_xe[mass][m.name] = xe.sigma_limit(mass=mass, 
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[xe.name])
        sigma_ge[mass][m.name] = ge.sigma_limit(mass=mass, 
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[ge.name])
        sigma_kimsiodLO[mass][m.name] = kimsiodLO.sigma_limit(mass=mass, 
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[kimsiodLO.name])
        sigma_kimsiodHI[mass][m.name] = kimsiodHI.sigma_limit(mass=mass, 
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[kimsiodHI.name])
        sigma_picoA[mass][m.name] = picoA.sigma_limit(mass=mass, 
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[picoA.name])
        sigma_picoB[mass][m.name] = picoB.sigma_limit(mass=mass, 
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[picoB.name])
        sigma_picoC[mass][m.name] = picoC.sigma_limit(mass=mass, 
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[picoC.name])
        sigma_picoD[mass][m.name] = picoD.sigma_limit(mass=mass, 
                                                    sigma_name=sigma_names[m.name],
                                                    fnfp_name=fnfp_names[m.name],
                                                    fnfp_val=fnfp_vals[m.name],
                                                    Nbackground=Nbg[picoD.name])
        sigma_vals[mass][m.name] = min(sigma_ge[mass][m.name],
                                     sigma_xe[mass][m.name],
                                     sigma_kimsiodLO[mass][m.name],
                                     sigma_kimsiodHI[mass][m.name],
                                     sigma_picoA[mass][m.name],
                                     sigma_picoB[mass][m.name],
                                     sigma_picoC[mass][m.name],
                                     sigma_picoD[mass][m.name])
        #print mass,m.name,sigma_vals[mass][m.name]

if DUMP_SIGMA_LIMS:
    fout = open(LIMFILE, 'wb')
    pickle.dump(sigma_vals, fout)
    
    
        
if DO_SIM:
    cmds = []
    count = 0
    for experiment in SINGLE_EXPERIMENTS:
        for mass in MASSES:
            for simmod in SIMMODELS:
                sigma_name = simmod.param_names[1]
                for i in range(STARTSIM,NSIM+STARTSIM):
                    simname='sim%d' % i
                    cmd = SCRIPTS_PATH + 'UVrunner.py --simname {} --simpars mass {} --parvals {} {:.16f} -e {}'.format(simname, sigma_name, mass, sigma_vals[mass][simmod.name], experiment)
                    if len(simmod.fixed_params) > 0:
                        cmd += ' --fixedsimnames {} --fixedsimvals {}'.format(simmod.fixed_params.keys()[0], simmod.fixed_params.values()[0])
                    cmds.append(cmd)
                    count += 1

    print  '\n There will be {} sims.\n'.format(count)
    if count < NSIMGROUPS:
        NSIMGROUPS = count

    for i in range(NSIMGROUPS):
        fout=open('runs_uv/simUVcommands_{}_{}.sh'.format(TAG, i+1), 'w')
        for cmd in cmds[i::NSIMGROUPS]:
            fout.write('{}\n'.format(cmd))
        fout.close()

    fout = open('runs_uv/simallUV_commandrunner_{}.sh'.format(TAG), 'w')
    fout.write('#! /bin/bash\n')
    fout.write('#$ -l h_rt=12:00:00\n')
    fout.write('#$ -cwd\n')
    fout.write('#$ -t 1-{}\n'.format(NSIMGROUPS))
    fout.write('#$ -V\n')
    fout.write('bash simUVcommands_{}_$SGE_TASK_ID.sh\n'.format(TAG))
    fout.close()
                                




if DO_FIT:
    cmds = []
    count = 0
    for experiment in EXPERIMENTS:
        for mass in MASSES:
            for simmod in SIMMODELS:
                sigma_name = simmod.param_names[1]

                for i in range(STARTFIT,NFIT+STARTFIT):
                    simname='sim%d' % i
                    cmd = SCRIPTS_PATH + 'UVrunner.py --fit --vis --simmodelname {} --simname {} --simpars mass {} --parvals {} {:.16f} -e {}'.format(simmod.name,
                                                                                                                                                 simname, 
                                                                                                                                                 sigma_name, 
                                                                                                                                                 mass, 
                                                                                                                                                 sigma_vals[mass][simmod.name], 
                                                                                                                                                 experiment)
                    if len(simmod.fixed_params) > 0:
                        cmd += ' --fixedsimnames {} --fixedsimvals {}'.format(simmod.fixed_params.keys()[0], simmod.fixed_params.values()[0])
                        
                    for fitmod in FITMODELS:
                        cmdfull = cmd + ' --fitmodelname {} --fitpars mass {}  --prior {}'.format(fitmod.name,fitmod.param_names[1],PRIOR_TYPE)
                        if len(fitmod.fixed_params) > 0:
                            cmdfull +=  ' --fixedfitnames {}  --fixedfitvals {}'.format(fitmod.fixed_params.keys()[0], fitmod.fixed_params.values()[0])
                        cmds.append(cmdfull)
                        count += 1
                                

    print  '\n There will be {} runs.\n'.format(count)
    if count < NGROUPS:
        NGROUPS = count

    for i in range(NGROUPS):
        fout=open('runs_uv/UVcommands_{}_{}.sh'.format(TAG, i+1), 'w')
        for cmd in cmds[i::NGROUPS]:
            fout.write('{}\n'.format(cmd))
        fout.close()

    fout = open('runs_uv/allUV_commandrunner_{}.sh'.format(TAG), 'w')
    fout.write('#! /bin/bash\n')
    fout.write('#$ -l h_rt=12:00:00\n')
    fout.write('#$ -cwd\n')
    fout.write('#$ -t 1-{}\n'.format(NGROUPS))
    fout.write('#$ -V\n')
    fout.write('bash UVcommands_{}_$SGE_TASK_ID.sh\n'.format(TAG))
    fout.close()

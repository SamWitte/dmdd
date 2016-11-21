#!/usr/bin/env python
#NOTE: this script expects that sims already exist. if some happen to not be there, it will create them; however, if all are absent, then there will be file-creation conflicts for sims the way this runner is set up.

import subprocess as sp
import os,sys,fnmatch
import argparse
import pickle

import dmdd
#import eff

parser = argparse.ArgumentParser()
parser.add_argument('--nofit', action='store_true')
parser.add_argument('--dosim', action='store_true')
parser.add_argument('--dumplims', action='store_true')
parser.add_argument('--limfile',default=None)
parser.add_argument('--tag',default='')
parser.add_argument('--masses',nargs='+',default=[20],type=float)
parser.add_argument('--ngroups', type=int, default=512)
parser.add_argument('--nsimgroups', type=int, default=512)
parser.add_argument('--nsim', type=int, default=1)
parser.add_argument('--startsim', type=int, default=1)
parser.add_argument('--nfit', type=int, default=1)
parser.add_argument('--startfit', type=int, default=1)
parser.add_argument('--prior',default='logflat')
parser.add_argument('--experiments',nargs='+',default=['Xe'])#['F','Ge', 'Xe','Ge Xe','I','Ge Xe I','Ge Xe F'],['Ilo','Xelo', 'Xehi','Xewide'],['Ar','Ge Xe Ar'],['He', 'Na', 'Ge','Ge He','Ge Na']
parser.add_argument('--path',default=os.environ['DMDD_AM_MAIN_PATH']+'/../')
parser.add_argument('--simtime',default=True)
parser.add_argument('--simGF', default=False)
parser.add_argument('--fittime',default=False)
parser.add_argument('--fitGF', default=False)
parser.add_argument('--analyze_energy', default=True)
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
SCRIPTS_PATH = '/home/verag/Projects/Repositories/dmdd-AM/'#args.path
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

SI_Higgs_sim = dmdd.UV_Model('SI_Higgs', ['mass', 'sigma_si'], fixed_params={'fnfp_si': 1}, time_info=args.simtime)
anapole_sim = dmdd.UV_Model('Anapole', ['mass','sigma_anapole'], time_info=args.simtime)
SI_Higgs_fit = dmdd.UV_Model('SI_Higgs', ['mass', 'sigma_si'], fixed_params={'fnfp_si': 1}, time_info=args.fittime)
anapole_fit = dmdd.UV_Model('Anapole', ['mass','sigma_anapole'], time_info=args.fittime)

SIMMODELS = [SI_Higgs_sim]
FITMODELS = [SI_Higgs_fit, anapole_fit]


##get upper limits for a given mass:
pandax=dmdd.Experiment('PandaX','xenon', 4., 30., 188., dmdd.eff.efficiency_Xe, 0., 1., energy_resolution=True)
lux=dmdd.Experiment('LUX','xenon', 5., 50., 70., dmdd.eff.efficiency_Xe, 0., 1., energy_resolution=True)
cdmslite=dmdd.Experiment('CDMSlite','germanium', 0.840, 6, 0.0164, dmdd.eff.efficiency_Xe, 0., 1., energy_resolution=True)
supercdms=dmdd.Experiment('SuperCDMS','germanium', 1.6, 12, 1.581, dmdd.eff.efficiency_Xe, 0., 1., energy_resolution=True)
sigma_vals = {}
sigma_pandax={}
sigma_lux={}
sigma_cdmslite={}
sigma_supercdms={}
for mass in MASSES:
    sigma_vals[mass] = {}
    sigma_lux[mass]={}
    sigma_cdmslite[mass]={}
    sigma_supercdms[mass]={}
    sigma_pandax[mass]={}
    for m in SIMMODELS:
        if len(m.fixed_params) > 0:
            sigma_lux[mass][m.name] = lux.sigma_limit(mass=mass, sigma_name=m.param_names[1],
                                                       fnfp_name=m.fixed_params.keys()[0],
                                                       fnfp_val=m.fixed_params.values()[0],
                                                       Nbackground=3.6)
            sigma_cdmslite[mass][m.name] = cdmslite.sigma_limit(mass=mass, sigma_name=m.param_names[1],
                                                       fnfp_name=m.fixed_params.keys()[0],
                                                       fnfp_val=m.fixed_params.values()[0])
            sigma_supercdms[mass][m.name] = supercdms.sigma_limit(mass=mass, sigma_name=m.param_names[1],
                                                       fnfp_name=m.fixed_params.keys()[0],
                                                       fnfp_val=m.fixed_params.values()[0])
            sigma_pandax[mass][m.name] = pandax.sigma_limit(mass=mass, sigma_name=m.param_names[1],
                                                            fnfp_name=m.fixed_params.keys()[0],
                                                            fnfp_val=m.fixed_params.values()[0],
                                                            Nbackground=13.95)

        else:
            sigma_lux[mass][m.name] = lux.sigma_limit(mass=mass, sigma_name=m.param_names[1], Nbackground=3.6)
            sigma_cdmslite[mass][m.name] = cdmslite.sigma_limit(mass=mass, sigma_name=m.param_names[1], Nbackground=4.)
            sigma_supercdms[mass][m.name] = supercdms.sigma_limit(mass=mass, sigma_name=m.param_names[1], Nbackground=4.)
            sigma_pandax[mass][m.name] = pandax.sigma_limit(mass=mass, sigma_name=m.param_names[1],
                                                            Nbackground=13.95)
        sigma_vals[mass][m.name]=min(sigma_supercdms[mass][m.name],
                                     sigma_cdmslite[mass][m.name],
                                     sigma_lux[mass][m.name],
                                     sigma_pandax[mass][m.name])
       
#####
if DUMP_SIGMA_LIMS:
    fout = open(LIMFILE, 'wb')
    pickle.dump(sigma_vals, fout)
    fout.close()

if DO_SIM:
    cmds = []
    count = 0
    for experiment in SINGLE_EXPERIMENTS:
        
        for mass in MASSES:
            for simmod in SIMMODELS:
                sigma_name = simmod.param_names[1]
                for i in range(STARTSIM,NSIM+STARTSIM):
                    simname='sim%d' % i
                    cmd = SCRIPTS_PATH+'UVrunner.py ' +\
                          '--simname {} --simpars mass {} '.format(simname, sigma_name) +\
                          ' --parvals {} {:.16f} '.format(mass,sigma_vals[mass][simmod.name]) +\
                          '-e {} --simGF {} --simtime {}'.format(experiment, args.simGF, args.simtime)
                          
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
                    cmd = SCRIPTS_PATH+'UVrunner.py --fit --vis --simmodelname {} '.format(simmod.name) +\
                            '--simname {} --simpars mass {} '.format(simname,sigma_name) +\
                            '--parvals {} {:.16f} '.format(mass, sigma_vals[mass][simmod.name]) +\
                            '-e {} --simtime {} --simGF {} --fittime {} --fitGF {} --analyze_energy {}'.format(experiment, args.simtime, args.simGF, args.fittime, args.fitGF, args.analyze_energy)
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

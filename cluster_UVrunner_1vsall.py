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
parser.add_argument('--masses',nargs='+',default=[500],type=float)
parser.add_argument('--ngroups', type=int, default=512)
parser.add_argument('--nsimgroups', type=int, default=512)
parser.add_argument('--nsim', type=int, default=1)
parser.add_argument('--startsim', type=int, default=1)
parser.add_argument('--nfit', type=int, default=1)
parser.add_argument('--startfit', type=int, default=1)
parser.add_argument('--prior',default='logflat')
parser.add_argument('--experiments',nargs='+',default=['Xe'])#['F','Ge', 'Xe','Ge Xe','I','Ge Xe I','Ge Xe F'],['Ilo','Xelo', 'Xehi','Xewide'],['Ar','Ge Xe Ar'],['He', 'Na', 'Ge','Ge He','Ge Na']
parser.add_argument('--path',default='/Users/SamWitte/Desktop/dmdd/')
parser.add_argument('--time',default=True)
#'/Users/verag/Research/Repositories/dmdd_2014/scripts/' #macbook


args = parser.parse_args()
DO_FIT = not args.nofit
DO_SIM = True
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

SI_Higgs = dmdd.UV_Model('SI_Higgs', ['mass', 'sigma_si'], fixed_params={'fnfp_si': 1}, time_info=args.time)
millicharge = dmdd.UV_Model('Millicharge', ['mass', 'sigma_si_massless'], fixed_params={'fnfp_si_massless': 0}, time_info=args.time)

SD_flavoruniversal = dmdd.UV_Model('SD_fu', ['mass','sigma_sd'], fixed_params={'fnfp_sd': -1.1}, time_info=args.time)
SD_Zmediated = dmdd.UV_Model('SD_Z-mediated', ['mass','sigma_sd'], fixed_params={'fnfp_sd': -13.2}, time_info=args.time)
SD_Moira = dmdd.UV_Model('SI_Moira', ['mass','sigma_sd'], fixed_params={'fnfp_sd': 0.}, time_info=args.time)

anapole = dmdd.UV_Model('Anapole', ['mass','sigma_anapole'], time_info=args.time)
magdip_heavy = dmdd.UV_Model('Mag.dip.heavy', ['mass','sigma_magdip'], time_info=args.time)
magdip_0 = dmdd.UV_Model('Mag.dip.light', ['mass','sigma_magdip_massless'], time_info=args.time)
elecdip_heavy = dmdd.UV_Model('Elec.dip.heavy', ['mass','sigma_elecdip'], time_info=args.time)
elecdip_0 = dmdd.UV_Model('Elec.dip.light', ['mass','sigma_elecdip_massless'], time_info=args.time)

f1 = dmdd.UV_Model('f1', ['mass','sigma_f1'], fixed_params={'fnfp_f1': 1.}, time_info=args.time)
f2_Higgs = dmdd.UV_Model('f2_Higgs', ['mass','sigma_f2'], fixed_params={'fnfp_f2': -0.05}, time_info=args.time)
f2_flavoruniversal = dmdd.UV_Model('f2_flavor-universal', ['mass','sigma_f2'], fixed_params={'fnfp_f2': 1.}, time_info=args.time)
f3_Higgs = dmdd.UV_Model('f3_Higgs', ['mass','sigma_f3'], fixed_params={'fnfp_f3': -0.05}, time_info=args.time)
f3_flavoruniversal = dmdd.UV_Model('f3_flavor-universal', ['mass','sigma_f3'], fixed_params={'fnfp_f3': 1.}, time_info=args.time)

LS = dmdd.UV_Model('LS', ['mass','sigma_LS'], fixed_params={'fnfp_LS': 0.}, time_info=args.time)



ALLMODELS = [SI_Higgs, millicharge, SD_flavoruniversal, SD_Zmediated, SD_Moira, anapole, magdip_heavy, magdip_0, elecdip_heavy, elecdip_0, f1, f2_Higgs, f2_flavoruniversal, f3_Higgs, f3_flavoruniversal]

MODELS1 = [SI_Higgs, anapole]

SIMMODELS = [SI_Higgs] #[SI_Higgs,elecdip_heavy,elecdip_0]
FITMODELS = MODELS1


##get upper limits for a given mass:
lux=dmdd.Experiment('LUX','xenon', 5, 23, 30.7, dmdd.eff.efficiency_Xe, 0., 1., energy_resolution=True)
cdmslite=dmdd.Experiment('CDMSlite','germanium', 0.840, 6, 0.0164, dmdd.eff.efficiency_Xe, 0., 1., energy_resolution=True)
supercdms=dmdd.Experiment('SuperCDMS','germanium', 1.6, 12, 1.581, dmdd.eff.efficiency_Xe, 0., 1., energy_resolution=True)
sigma_vals = {}
sigma_lux={}
sigma_cdmslite={}
sigma_supercdms={}
for mass in MASSES:
    sigma_vals[mass] = {}
    sigma_lux[mass]={}
    sigma_cdmslite[mass]={}
    sigma_supercdms[mass]={}
    for m in SIMMODELS:
        if len(m.fixed_params) > 0:
            sigma_lux[mass][m.name] = lux.sigma_limit(mass=mass, sigma_name=m.param_names[1],
                                                       fnfp_name=m.fixed_params.keys()[0],
                                                       fnfp_val=m.fixed_params.values()[0])
            sigma_cdmslite[mass][m.name] = cdmslite.sigma_limit(mass=mass, sigma_name=m.param_names[1],
                                                       fnfp_name=m.fixed_params.keys()[0],
                                                       fnfp_val=m.fixed_params.values()[0])
            sigma_supercdms[mass][m.name] = supercdms.sigma_limit(mass=mass, sigma_name=m.param_names[1],
                                                       fnfp_name=m.fixed_params.keys()[0],
                                                       fnfp_val=m.fixed_params.values()[0])

        else:
            sigma_lux[mass][m.name] = lux.sigma_limit(mass=mass, sigma_name=m.param_names[1])
            sigma_cdmslite[mass][m.name] = cdmslite.sigma_limit(mass=mass, sigma_name=m.param_names[1])
            sigma_supercdms[mass][m.name] = supercdms.sigma_limit(mass=mass, sigma_name=m.param_names[1])
        sigma_vals[mass][m.name]=min(sigma_supercdms[mass][m.name],
                                     sigma_cdmslite[mass][m.name],
                                     sigma_lux[mass][m.name])
       
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
                    cmd = 'cd '+ SCRIPTS_PATH + '\n' + 'python UVrunner.py --simname {} --simpars mass {} --parvals {} {:.16f} -e {}'.format(simname, sigma_name, mass, 
                                                                                                                                                       sigma_vals[mass][simmod.name],
                                                                                                                                                       experiment)
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
                    cmd = 'cd ' + SCRIPTS_PATH + '\n' + 'python UVrunner.py --fit --vis --simmodelname {} --simname {} --simpars mass {} --parvals {} {:.16f} -e {} --time {}'.format(simmod.name,
                                                                                                                                                 simname, 
                                                                                                                                                 sigma_name, 
                                                                                                                                                 mass, 
                                                                                                                                                 sigma_vals[mass][simmod.name], 
                                                                                                                                                 experiment, args.time)
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

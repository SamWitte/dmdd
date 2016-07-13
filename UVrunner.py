#!/usr/bin/env python

"""This is setup only to run rate_Haxton models.
usage example (from scripts directory):
./UVrunner.py --simname test --simpars mass sigma_si --parvals 50. 10. -e Ge Xe --fitpars sigma_sd mass --prior logflat --vis
"""
import time
#start = time.time()

import matplotlib 
matplotlib.use('agg')
import argparse


import dmdd
#import eff
#import rate_Haxton
from experiment_parameters import *
from multinest_parameters import *


parser = argparse.ArgumentParser()

parser.add_argument('--simmodelname',default='simmodel')
parser.add_argument('--fitmodelname',default='fitmodel')
parser.add_argument('--simpars',nargs='+',default=['mass'])
parser.add_argument('--fitpars',nargs='+',default=['mass'])
parser.add_argument('--parvals',nargs='+',default=[50.],type=float)

parser.add_argument('--fixedsimnames',nargs='+',default=[])
parser.add_argument('--fixedfitnames',nargs='+',default=[])
parser.add_argument('--fixedsimvals',nargs='+',default=[],type=float)
parser.add_argument('--fixedfitvals',nargs='+',default=[],type=float)


parser.add_argument('--simname',default='sim1')
parser.add_argument('-e','--exps',nargs='+',default=['Xe'])

parser.add_argument('--prior',default='logflat')

parser.add_argument('--asimov', action='store_true')
parser.add_argument('--nasbin', type=int, default=20)
parser.add_argument('--forcesim', action='store_true')
parser.add_argument('--vis', action='store_true')
parser.add_argument('--fit', action='store_true')

parser.add_argument('--nlive', type=int, default=n_live_points)
parser.add_argument('--evtol',type=float, default=evidence_tolerance)
parser.add_argument('--seff',type=float, default=sampling_efficiency)
parser.add_argument('--resume', action='store_true')
parser.add_argument('--base', default=basename)
parser.add_argument('--time', default='T')
parser.add_argument('--GF', default='F')

args = parser.parse_args()

if args.GF == 'T':
    GF = True
elif args.GF == 'F':
    GF = False

experiments = []
for i,experiment in enumerate(args.exps):
    efficiency = dmdd.eff.efficiency_Xe
    experiments.append(dmdd.Experiment(experiment,Target[experiment],
                                     Qmin[experiment], Qmax[experiment],
                                     Exposure[experiment], efficiency, 0., 1.,
                                     energy_resolution=EnergyResolution[experiment]))



fixedsim_params = {}
for i,par in enumerate(args.fixedsimnames):
    fixedsim_params[par] = args.fixedsimvals[i]

fixedfit_params = {}
for i,par in enumerate(args.fixedfitnames):
    fixedfit_params[par] = args.fixedfitvals[i]

simmodel = dmdd.UV_Model(args.simmodelname, args.simpars, fixed_params=fixedsim_params, time_info=args.time, GF=GF)
fitmodel = dmdd.UV_Model(args.fitmodelname, args.fitpars, fixed_params=fixedfit_params, time_info=args.time, GF=GF)

param_values = {}
for i,par in enumerate(args.simpars):
    param_values[par] = args.parvals[i]

print 'in UVrunner: ', param_values
Haxton_run = dmdd.MultinestRun(args.simname, experiments, simmodel, param_values, fitmodel,
                                force_sim=args.forcesim, prior_ranges=prior_ranges, prior=args.prior,
                                asimov=args.asimov, nbins_asimov=args.nasbin,
                                n_live_points=args.nlive, evidence_tolerance=args.evtol,
                                sampling_efficiency=args.seff, resume=args.resume, basename=args.base,
                                time_info=args.time, GF=GF, TIMEONLY=False)
if args.fit:
    Haxton_run.fit()
if args.vis:
    Haxton_run.visualize()

"""
end = time.time()
print '\n timing: {:.2f}min (sim: {}, {} fit:{}, {} mass: {})\n'.format((end - start) / 60.,
                                                                        args.simpars[1],
                                                                        args.fixedsimnames[0],
                                                                        args.fitpars[1],
                                                                        args.fixedfitnames[0],
                                                                        args.simpars[0])
"""

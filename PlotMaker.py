import UVplots
from UVplots import *

line_plots_1vsall(nsim=1, startsim=1, masses=[50.],
                      experiment_names=['Ge','Xe'],
                      simmodels=[anapole], models=[SI_Higgs, anapole], time_info='Both',
                      GF=True, filelabel='TEST', allverbose=True,
                      verbose=True,
                      results_root='/Users/SamWitte/Desktop/dmdd/Storage/results_uv/',
                      saveplots=True, alpha=0.3, xoffset=0.1, fs=20, fs2=18, sigma_lim_file=None,
                      colors_list=['Aqua','Red'])


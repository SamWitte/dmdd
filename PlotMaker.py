import UVplots
from UVplots import *

line_plots_1vsall(nsim=50, startsim=1, masses=[500.],
                      experiment_names=['Xe'],
                      simmodels=[SI_Higgs], models=[SI_Higgs, anapole], time_info='Both',
                      filelabel='Xe_Time_vs_NoTime', allverbose=True,
                      verbose=True,
                      results_root='/Users/SamWitte/Desktop/dmdd/Storage/results_uv/',
                      saveplots=True, alpha=0.3, xoffset=0.1, fs=20, fs2=18, sigma_lim_file=None,
                      colors_list=['Aqua','Red'])


import UVplots
from UVplots import *

OneDhistogram(nsim=50, startsim=1, masses=[500.],
                      experiment_names=['Xe', 'XeTriple', 'Xe10x'],
                      simmodels=[anapole], models=[SI_Higgs, anapole], time_info='Both',
                      GF=False, hspace = (1.06 * 50. ** (-1. / 5.)), filelabel='TEST_PDF', allverbose=True,
                      verbose=True,
                      results_root=os.environ['DMDD_AM_MAIN_PATH']+'/results_uv_noGF/',
                      saveplots=True, alpha=0.2, xoffset=0.1, fs=20, fs2=18, sigma_lim_file=None,
                      colors_list=['Aqua','Red','Black','Magenta','Green','Orange'])

OneDhistogram_timeDiff(nsim=50, startsim=1, masses=[500.],
                  experiment_names=['Xe', 'XeTriple', 'Xe10x'],
                  simmodels=[anapole], models=[SI_Higgs, anapole], GF=False,
                  hspace = (1.06 * 50. ** (-1. / 5.)), filelabel='TEST', allverbose=True, verbose=True,
                  results_root=os.environ['DMDD_AM_MAIN_PATH']+'/results_uv_noGF/', timeonly=False,
                  saveplots=True, alpha=0.3, xoffset=0.1, fs=16, fs2=18, sigma_lim_file=None,
                  colors_list=['Aqua','Red','Black','Green','Magenta','Orange'])

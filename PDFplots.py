import UVplots
from UVplots import *

OneDhistogram(nsim=50, startsim=1, masses=[20.],
                      experiment_names=['Ge+Xe'],
                      simmodels=[anapole], models=[SI_Higgs, anapole], time_info='Both',
                      GF=True, hspace = .15, filelabel='TEST_PDF', allverbose=True,
                      verbose=True,
                      results_root=os.environ['DMDD_AM_MAIN_PATH']+'/results_uv/',
                      saveplots=True, alpha=0.2, xoffset=0.1, fs=20, fs2=18, sigma_lim_file=None,
                      colors_list=['Aqua','Red','Black','Magenta','Green','Orange'])


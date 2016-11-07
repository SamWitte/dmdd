import numpy as np
import os
cimport numpy as np
cimport cython
from cpython cimport bool
#from scipy.interpolate import griddata,interp1d,interp2d


DTYPE = np.float
ctypedef np.float_t DTYPE_t
cdef DTYPE_t pi = np.pi #3.14159265359

eta0_a0_tabbed = np.loadtxt(os.environ['DMDD_AM_MAIN_PATH']+'/..'+'/dmdd/eta0_a0.dat')
#cdef np.float_t[:,:] eta0_a0_tabbed = np.loadtxt(os.environ['DMDD_AM_MAIN_PATH']+'/..'+'/dmdd/eta0_a0.dat')
eta0_a1_tabbed = np.loadtxt(os.environ['DMDD_AM_MAIN_PATH']+'/..'+'/dmdd/eta0_a1.dat')
#cdef np.float_t[:,:] eta0_a1_tabbed = np.loadtxt(os.environ['DMDD_AM_MAIN_PATH']+'/..'+'/dmdd/eta0_a1.dat')
eta0_b1_tabbed = np.loadtxt(os.environ['DMDD_AM_MAIN_PATH']+'/..'+'/dmdd/eta0_b1.dat')
#cdef np.float_t[:,:] eta0_b1_tabbed = np.loadtxt(os.environ['DMDD_AM_MAIN_PATH']+'/..'+'/dmdd/eta0_b1.dat')
eta1_a0_tabbed = np.loadtxt(os.environ['DMDD_AM_MAIN_PATH']+'/../'+'/dmdd/eta1_a0.dat')
#cdef np.float_t[:,:] eta1_a0_tabbed = np.loadtxt(os.environ['DMDD_AM_MAIN_PATH']+'/../'+'/dmdd/eta1_a0.dat')
eta1_a1_tabbed = np.loadtxt(os.environ['DMDD_AM_MAIN_PATH']+'/../'+'/dmdd/eta1_a1.dat')
#cdef np.float_t[:,:] eta1_a1_tabbed = np.loadtxt(os.environ['DMDD_AM_MAIN_PATH']+'/../'+'/dmdd/eta1_a1.dat')
eta1_b1_tabbed = np.loadtxt(os.environ['DMDD_AM_MAIN_PATH']+'/../'+'/dmdd/eta1_b1.dat')
#cdef np.float_t[:,:] eta1_b1_tabbed = np.loadtxt(os.environ['DMDD_AM_MAIN_PATH']+'/../'+'/dmdd/eta1_b1.dat')
 

cdef extern from "math.h":
    double sin(double)
    double cos(double)
    double tanh(double)
    double sqrt(double)
    double atan2(double,double)
    double acos(double)
    double abs(double)
    double log(double)
    double ceil(double)
    double fabs(double)
    double exp(double)
    double erf(double)
    double log10(double)


def trapz(np.ndarray[DTYPE_t] y, np.ndarray[DTYPE_t] x):
    """
    Integrator function, the same as numpy.trapz.
    The points y,x must be given in increasing order.
    """
    if len(x) != len(y):
        raise ValueError('x and y must be same length')
    cdef long npts = len(x)
    cdef DTYPE_t tot = 0
    cdef unsigned int i
    
    for i in range(npts-1):
        tot += 0.5*(y[i]+y[i+1])*(x[i+1]-x[i])
    return tot



def interp1d(np.ndarray[DTYPE_t] x, np.ndarray[DTYPE_t] y, DTYPE_t x0):
    cdef int i = 0
    cdef DTYPE_t res

    while x[i] < x0:
        i = i + 1
        
    if x[i] == x0:
        return y[i]
    if i == 0:
        return y[i] - (y[i+1] - y[i]) * (x0 - x[i]) / (x[i] - x[i])
        
    res = y[i-1] + (y[i] - y[i-1]) * (x0 - x[i-1]) / (x[i] - x[i-1])
    return res



def log_fact(DTYPE_t xx):
    """
    Approximate formula for ln(x!).
    """
    cdef DTYPE_t ser=1.000000000190015
    cdef DTYPE_t x,y,tmp
    cof = [76.18009172947146,-86.50532032941677,
           24.01409824083091,-1.231739572450155,
           0.1208650973866179e-2,-0.5395239384953e-5]
    
    y=xx+1
    x=xx+1
    tmp=x+5.5
    tmp -= (x+0.5)*log(tmp)
    for c in cof:
        y=y+1
        ser += c/y
    return -tmp+log(2.5066282746310005*ser/x)


def eta(DTYPE_t v_min, DTYPE_t v_esc, DTYPE_t v_rms, DTYPE_t v_lag):
    """
    This is the correctly scaled velocity integral for a rate
    with no special velocity dependence.

    The input units for all velocities are km/s.
    """
    cdef DTYPE_t res,xmin,xe,xesc, integral,norm
    xmin = v_min/v_rms
    xe = v_lag/v_rms + 17.72/v_rms
    xesc = v_esc/v_rms
    norm = 1./(4*xe*(sqrt(pi)*0.5*erf(xesc)-exp(-xesc**2)*xesc*(1+2*xesc**2/3.)))

    if xesc - xe - xmin > 0:
        integral = (sqrt(pi)*(erf(xe+xmin)-erf(xmin-xe))-4*xe*exp(-xesc**2)*(1+xesc**2-xmin**2-xe**2./3))*norm
    elif xesc - xmin + xe > 0 and -xesc+xmin+xe>0:
        integral =  (sqrt(pi)*(erf(xesc)-erf(xmin-xe))-2*exp(-xesc**2)*(xesc+xe-xmin-1./3*(xe-2*xesc-xmin)*(xesc+xe-xmin)**2))*norm
    else:
        integral =  0
    
    res = 3*10**5*integral/v_rms
    
    return res


def eta_GF(DTYPE_t v_min, DTYPE_t time, bool time_info):
    """
    Calculated only for v_esc = 533km/s, v_rms = v_lag = 220 km/s
    
    Also, if vmin_
    """
    cdef DTYPE_t res
    cdef DTYPE_t vmin_max=700.
    cdef DTYPE_t vlag = 220.
    cdef DTYPE_t mod_amp = 29.8 * 0.49
    cdef DTYPE_t mod_phase = 0.42
    cdef DTYPE_t twopi = 2. * pi


    global eta0_a0_tabbed
    global eta0_a1_tabbed
    global eta0_b1_tabbed
    
    cdef np.ndarray[DTYPE_t] a0_x = eta0_a0_tabbed[:,0]
    cdef np.ndarray[DTYPE_t] a0_y = eta0_a0_tabbed[:,1]
    cdef np.ndarray[DTYPE_t] a1_x = eta0_a1_tabbed[:,0]
    cdef np.ndarray[DTYPE_t] a1_y = eta0_a1_tabbed[:,1]
    cdef np.ndarray[DTYPE_t] b1_x = eta0_b1_tabbed[:,0]
    cdef np.ndarray[DTYPE_t] b1_y = eta0_b1_tabbed[:,1]
    

    if v_min <= vmin_max:
        if time_info:
            res = 3. * 10.**5. * (interp1d(a0_x, a0_y, v_min) + 
                  interp1d(a1_x, a1_y, v_min) * cos(2. * pi * (time - 0.4178)) +
                  interp1d(b1_x, b1_y, v_min) * sin(2. * pi * (time - 0.4178))) 
        else:   
            res = 3. * 10.**5. * (interp1d(a0_x, a0_y, v_min))
    else:
        if time_info:
            vlag += mod_amp * cos(twopi * (time - mod_phase))
        res = eta(v_min, 533., 220., vlag)
    
    return res

def zeta_GF(DTYPE_t v_min, DTYPE_t time, bool time_info):
    """
    This is the correctly scaled velocity integral for a rate
    with no special velocity dependence.

    The input units for all velocities are km/s.
    """
    cdef DTYPE_t res
    cdef DTYPE_t vmin_max=700.
    cdef DTYPE_t vlag = 220.
    cdef DTYPE_t mod_amp = 29.8 * 0.49
    cdef DTYPE_t mod_phase = 0.42
    cdef DTYPE_t twopi = 2. * pi

    global eta1_a0_tabbed
    global eta1_a1_tabbed
    global eta1_b1_tabbed
    
    cdef np.ndarray[DTYPE_t] a0_x = eta1_a0_tabbed[:,0]
    cdef np.ndarray[DTYPE_t] a0_y = eta1_a0_tabbed[:,1]
    cdef np.ndarray[DTYPE_t] a1_x = eta1_a1_tabbed[:,0]
    cdef np.ndarray[DTYPE_t] a1_y = eta1_a1_tabbed[:,1]
    cdef np.ndarray[DTYPE_t] b1_x = eta1_b1_tabbed[:,0]
    cdef np.ndarray[DTYPE_t] b1_y = eta1_b1_tabbed[:,1]
    

    if v_min <= vmin_max:
        if time_info:
            res = (interp1d(a0_x, a0_y, v_min) + 
                   interp1d(a1_x, a1_y, v_min) * cos(2. * pi * (time - 0.4178)) +
                   interp1d(b1_x, b1_y, v_min) * sin(2. * pi * (time - 0.4178)))  / (3.*10**5.) #+
        else:
            res = interp1d(a0_x, a0_y, v_min) / (3.*10**5.) 
    else:
        if time_info:
            vlag += mod_amp * cos(twopi * (time - mod_phase))
        res = zeta(v_min, 533., 220., vlag)

    return res


def zeta(DTYPE_t v_min,DTYPE_t v_esc,DTYPE_t v_rms, DTYPE_t v_lag):
    """
    This is the correctly scaled velocity integral for a rate with additional velocity^2 dependence. The input units for all velocities are km/s.
    """
    cdef DTYPE_t res,xmin,xe,xesc, integral,norm
    xmin = v_min/v_rms
    xe = v_lag/v_rms + 17.72/v_rms
    xesc = v_esc/v_rms
    norm = 1./(4*xe*(sqrt(pi)*0.5*erf(xesc)-exp(-xesc**2)*xesc*(1+2*xesc**2/3.)))
    
    if xesc-xmin-xe > 0:
        integral = (1./2*sqrt(pi)*(2*xe**2+1)*(erf(xmin+xe)-erf(xmin-xe)) + exp(-(xmin-xe)**2)*(xe+xmin) + exp(-(xmin+xe)**2)*(xe-xmin) - exp(-xesc**2)*(4*(1+xesc**2+xe**2./3)*xe-2*xe**5/15.+4*xe**3*xesc**2/3.-2*xe*(xmin**4-xesc**4)))*norm
    elif xesc-xmin+xe > 0 and -xesc+xmin+xe > 0:
        integral = 2*exp(-xesc**2)*(1./3*(xmin**3-(xesc+xe)**3)-1./4*exp(-xmin**2-xe**2+xesc**2)*(pi**.5*(2*xe**2+1)*exp(xmin**2+xe**2)*(erf(xmin-xe)-erf(xesc))-2*exp(2*xmin*xe)*(xmin+xe)+2*(xesc+2*xe)*exp(xmin**2+xe**2-xesc**2))-1./30*(-xe**5+10*xe**3*xesc**2+10*xe**2*(2*xesc**3+xmin**3)+15*xe*(xesc**4-xmin**4)+4*xesc**5-10*xesc**2*xmin**3+6*xmin**5))*norm
    else:
        integral =  0
    res = integral*v_rms/(3.*10**5.)
    return res




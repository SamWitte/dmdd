#!/usr/local/bin/python
import numpy as np
import time
import dmdd
def timetest(witht=True):
    qs = np.linspace(2,10,10)
    tim = np.linspace(0,1.,10)
    if witht:
        n = 1
    else:
        n = 100
    if witht:
        startt = time.time()
        for i in xrange(n):
            res = dmdd.rate_UV.loglikelihood_time(qs,0.,1.,tim,dmdd.eff.efficiency_unit,sigma_si=100)
        endt = time.time()
    else:
        startt = time.time()
        for i in xrange(n):
            res = dmdd.rate_UV.loglikelihood(qs,dmdd.eff.efficiency_unit,sigma_si=100)
        endt = time.time()
    print('{}\n'.format((endt-startt)/n))

timetest()
timetest(witht=False)

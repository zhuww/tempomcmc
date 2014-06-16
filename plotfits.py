from datatools.fitstools import fitsfile
import sys
from pylab import *
import numpy as np
from decimal import Decimal
from Coordinate import RA, Dec
secperday = 24*3600
pi = 3.141592653589793
#Msun = 1.98892e30
Tsun = 4.925490947e-6


class MCMCresult(object):
    def __init__(self, mcmcfile):
        self.fitsfile = fitsfile(mcmcfile, ignore_missing_end=True)
        #ff = fitsfile(sys.argv[1]+".mcmc", ignore_missing_end=True)
        SPCPAR = self.fitsfile.hdread('SPCPAR').split('|')
        self.SPCPAR = SPCPAR
        cols = self.fitsfile.colinfo()
        names = cols.names
        tb = self.fitsfile.gettable()[:-1]
        self.table = tb
        #ntb = ff.gettable(cols=sys.argv[2:])
        print tb.shape

        try:
            if 'PB' in SPCPAR:
                Pb = float(self.fitsfile.hdread('PB'))
            else:
                Pb = tb.field('PB')
            self.Pb = Pb
        except:pass
        try:
            SINI = tb.field('SINI')
            self.SINI = SINI
        except:pass
        try:
            a = tb.field('A1')
            self.A1 = a
        except:pass
        try:
            M2 = tb.field('M2')
            M1 = (Pb*secperday/2/pi*sqrt(Tsun*(M2*SINI)**3/a**3)-M2)
            self.M1 = M1
            self.M2 = M2
        except:pass
        chisq = tb.field('chisq')
        bestidx = np.argmin(chisq)
        self.chisq = chisq
        self.bestidx = bestidx

if __name__ == '__main__':
    ff = MCMCresult(sys.argv[1]+".mcmc")
    results = []
    bests = []
    labels = []
    for par in sys.argv[2:]:
        if par == 'M1':
            results.append(ff.M1)
        else:
            results.append(ff.table.field(par))

        if par in ff.SPCPAR:
            labels.append('%s - %s' % (par,ff.fitsfile.hdread(par))) 
        elif par in ['RAJ', 'DECJ']:
            labels.append('%s - %s (sec)' % (par,ff.fitsfile.hdread(par))) 
        else:
            labels.append(par)
if len(results) == 2:
    xlabel(labels[0])
    ylabel(labels[1])
    plot(results[0], results[1], '.')
    plot(results[0][ff.bestidx], results[1][ff.bestidx], 'ro')
    show()
elif len(results) == 1:
    xlabel(labels[0])
    ylabel('probability')
    try:
        n, bins, patches = hist(results[0], 50, normed=1, facecolor = 'blue', )
    except:
        print results[0]
    show()


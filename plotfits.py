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

ff = fitsfile(sys.argv[1]+".mcmc")

SPCPAR = ff.hdread('SPCPAR').split('|')
cols = ff.colinfo()
names = cols.names
tb = ff.gettable()[:-1]
#ntb = ff.gettable(cols=sys.argv[2:])
print tb.shape

try:
    if 'PB' in SPCPAR:
        Pb = float(ff.hdread('PB'))
    else:
        Pb = tb.field('PB')
except:pass
try:
    SINI = tb.field('SINI')
except:pass
try:
    a = tb.field('A1')
except:pass
try:
    M2 = tb.field('M2')
    M1 = (Pb*secperday/2/pi*sqrt(Tsun*(M2*SINI)**3/a**3)-M2)
except:pass
chisq = tb.field('chisq')
bestidx = np.argmin(chisq)

if __name__ == '__main__':
    results = []
    bests = []
    labels = []
    for par in sys.argv[2:]:
        if par == 'M1':
            results.append(M1)
        else:
            results.append(tb.field(par))

        if par in SPCPAR:
            labels.append('%s - %s' % (par,ff.hdread(par))) 
        elif par in ['RAJ', 'DECJ']:
            labels.append('%s - %s (sec)' % (par,ff.hdread(par))) 
        else:
            labels.append(par)
if len(results) == 2:
    xlabel(labels[0])
    ylabel(labels[1])
    plot(results[0], results[1], '.')
    plot(results[0][bestidx], results[1][bestidx], 'ro')
    show()
elif len(results) == 1:
    xlabel(labels[0])
    ylabel('probability')
    try:
        n, bins, patches = hist(results[0], 50, normed=1, facecolor = 'blue', )
    except:
        print results[0]
    show()


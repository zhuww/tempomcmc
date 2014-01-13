from pylab import *
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import cPickle as pickle
import sys
from Coordinate import RA, Dec
secperday = 24*3600


dict = pickle.load(open('bestpar.p', 'rb'))
best = dict['BEST']
plist = dict['parameters'] 
MChain = pickle.load(open('MChain.p','rb'))
MarkovChain = MChain['Chain']
print len(MarkovChain)
pi = 3.141592653589793
G = 6.673e-11
Msun = 1.98892e30
c = 2.99792458e8
twopi = 6.283185307179586
fac = 1.536e-16 
im2 = plist.index('M2')
ipb = plist.index('PB')
isini = plist.index('SINI')
ia = plist.index('A1')
ichisq = plist.index('chisq')
M2 = np.array([float(p[im2])*Msun for p in MarkovChain])
Pb = np.array([float(p[ipb])*secperday for p in MarkovChain])
SINI = np.array([float(p[isini]) for p in MarkovChain])
a = np.array([float(p[ia])*c for p in MarkovChain])
M1 = (Pb/2/pi*sqrt(G*(M2*SINI)**3/a**3)-M2)/Msun
M2 = M2/Msun
chisq = [p[ichisq] for p in MarkovChain]
bestidx = chisq.index(min(chisq))

if __name__ == '__main__':
#def plotpar(par1, par2)
    #par1 = sys.argv[1]
    #par2 = sys.argv[2]
    results = []
    bests = []
    for par in sys.argv[1:]:
        if par == 'M1':
            results.append(M1)
        elif par == 'M2':
            results.append(M2)
        elif par == 'SINI':
            results.append(SINI)
        elif par == 'RAJ':
            i = plist.index(par)
            results.append(np.array([RA(x[i]).in_unit_degree for x in MarkovChain]))
        elif par == 'DECJ':
            i = plist.index(par)
            results.append(np.array([Dec(x[i]).in_unit_degree for x in MarkovChain]))
        else:
            i = plist.index(par)
            results.append(np.array([float(x[i]) for x in MarkovChain]))



if len(results) == 2:
    xlabel(sys.argv[1])
    ylabel(sys.argv[2])
    plot(results[0], results[1], '.')
    plot(results[0][bestidx], results[1][bestidx], 'ro')
    show()
elif len(results) == 1:
    xlabel(sys.argv[1])
    ylabel('probability')
    try:
        n, bins, patches = hist(results[0], 50, normed=1, facecolor = 'blue', )
    except:
        print results[0]
    show()


from pylab import *
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pickle 
import sys
from Coordinate import RA, Dec



if __name__ == '__main__':
    par1 = sys.argv[1]
    par2 = sys.argv[2]
    dict = pickle.load(open('bestpar.p', 'r'))
    best = dict['BEST']
    plist = dict['parameters']
    MChain = pickle.load(open('MChain.p','r'))
    MarkovChain = MChain['Chain']
    print len(MarkovChain)
    i = plist.index(par1)
    j = plist.index(par2)
    best1 = [best[i]]
    best2 = [best[j]]
    if par1 == 'RAJ':
        X = [RA(x[i]).in_unit_degree for x in MarkovChain]
        best1 = [RA(best1[0]).in_unit_degree]
    elif par1 == 'DECJ':
        X = [Dec(x[i]).in_unit_degree for x in MarkovChain]
        best1 = [Dec(best1[0]).in_unit_degree]
    else:
        X = [float(x[i]) for x in MarkovChain]
    if par2 == 'RAJ':
        Y = [RA(x[j]).in_unit_degree for x in MarkovChain]
        best2 = [RA(best2[0]).in_unit_degree]
    elif par2 == 'DECJ':
        Y = [Dec(x[j]).in_unit_degree for x in MarkovChain]
        best2 = [Dec(best2[0]).in_unit_degree]
    else:
        Y = [float(x[j]) for x in MarkovChain]

    xlabel(par1)
    ylabel(par2)
    plot(X, Y, '.')
    plot(best1, best2, 'ro')
    show()


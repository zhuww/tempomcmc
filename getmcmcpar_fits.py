from tempo import *
import cPickle as pickle
import sys
from pylab import *
import numpy as np
from decimal import Decimal
from Coordinate import RA, Dec
import fitsio


secperday = 24*3600
pi = 3.141592653589793
#Msun = 1.98892e30
Tsun = 4.925490947e-6
parfile = sys.argv[1]

pf = PARfile(parfile)
try:
    PSRname = pf.PSRJ
except:
    PSRname = pf.PSR
ff = fitsio.FITS(PSRname + '.mcmc')
data = ff[-1].read()
hdr = ff[-1].read_header()
plist = ff[-1].get_colnames()
LongPar = hdr['SPCPAR'].split('|')
for p in plist:
    vals = data[p]
    val = np.mean(vals)
    err = np.std(vals)
    if p in ['RAJ', 'DECJ']:
        D, H, S = hdr[p].split(':')
        newS = float(S)+val
        if newS >= 0 and newS < 60:
            pf.__dict__[p][0] = ':'.join([D, H, str(newS)])
        elif newS < 0:
            pf.__dict__[p][0] = ':'.join([D, str(int(H)-1), str(newS+60.)])
        else:
            pf.__dict__[p][0] = ':'.join([D, str(int(H)+1), str(newS-60.)])
        pf.__dict__[p][1] = Decimal(str(err))
    elif p in LongPar:
        pf.__dict__[p][0] = Decimal(hdr[p]) + Decimal(str(val))
        pf.__dict__[p][1] = Decimal(str(err))
    else:
        pf.__dict__[p] = [Decimal(str(val)), Decimal(str(err))]


if 'PAASCNODE' in plist:
    pf.__dict__['KOM'] = pf.__dict__['PAASCNODE']
    pf.manifest.append('KOM')
    pf.parameters['PAASCNODE'] = 1
    pf.parameters['KOM'] = 1

if 'M2' in plist and 'SINI' in plist:
    if 'PB' in LongPar:
        Pb = float(hdr['PB'])
    else:
        Pb = data['PB'] 
    SINI = data['SINI']
    a = data['A1']
    M2 = data['M2']
    M1 = Pb*secperday/2/np.pi*np.sqrt(Tsun*(M2*SINI)**3/a**3)-M2
    KIN = np.arcsin(SINI)*180./np.pi
    for p in ['M1', 'KIN',]:
        vals = locals()[p]
        val = np.mean(vals)
        err = np.std(vals)
        pf.__dict__[p] = [Decimal(str(val)), Decimal(str(err))]
        pf.manifest.append(p)
        pf.parameters[p] = 1

pf.write('mcmcfitsresult.par')

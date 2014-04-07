from tempo import *
from TexTable import deluxetable
from round import TexStyle as SF
import numpy as np
from ordereddict import OrderedDict
import decimal
from tools.Coordinate import RA, Dec


secperyear = 3600*24*365.24218967
secperday = 3600 * 24
PI = np.pi

def M1(pf):
    G = float(6.673e-11)
    Msun = float(1.98892e30)
    Tsun = float(4.925490947)
    c = float(2.99792458e8)
    m2 = float(pf.M2[0])
    I = float(pf.KIN[0])/180*np.pi
    Pb = float(pf.PB[0])*secperday
    a = float(pf.A1[0])
    #result = sqrt(930.998*m2**3*Pb**2/a**3) - m2
    return (Pb/2/PI*(sqrt(Tsun*(m2*sin(I))**3/a**3))-m2)

def B(pf):
    return 3.2e19*np.sqrt(np.abs(float(pf.F1[0]/pf.F0[0]**3)))


def Age(pf):
    return np.abs(float(pf.F0[0]/pf.F1[0]/secperyear/2))

def aveDM(pf):
    DMX, DMXErr, DMXR1, DMXR2 = pf.dmxlist
    dmx = [DMX[x] for x in DMX]
    dmxerr = np.std([float(x) for x in dmx])
    return sum(dmx)/len(dmx) , Decimal(str(dmxerr))

def parseerror(val, err):
    s, d ,e = err.as_tuple()
    if int(d[0]) == 1:# and int(d[1]) >= 5:
    #if int(d[0]) == 1 and not int(d[1]) == 0:
    #if int(d[0]) == 1: #and int(d[1]) >= 5:
        errdig = Decimal((s,d[:2],e+len(d)-2))
        errstr = Decimal((0, d[:2], 0))
        val = str(val.quantize(errdig))
    else:
        errdig = Decimal((s,d[:1],e+len(d)-1))
        errstr = Decimal((0, d[:1], 0))
        val = str(val.quantize(errdig))
    if not val.find('E') == -1:
        v,e = val.split('E')
        result = r'$%s(%s)\times10^{%s}$' % (v,errstr, e)
    else:
        result = r'%s(%s)' % (val,errstr)
    return result


#m = model('1713.DMX.noOMDOT.par.t2')
m = model('mcmcresult.par')
#m = model('bestmcmc.t2.par')
#m1 = model('bestmcmc.par')
#m.XDOT = m1.XDOT
#print aveDM(m)
#sys.exit(0)

Observable = OrderedDict([
        ('RAJ',r'Right Ascension, $\alpha$ (J2000)'),
        ('DECJ',r'Declination, $\delta$ (J2000)'),
        ('F0',r'Spin Frequecy $\nu$~(s$^{-1}$)'),
        ('F1',r"Spin down rate $\nu'$ (s$^{-2}$)"),
        ('F2',r"Second spin frequency derivative $\nu''$ (s$^{-3}$)"),
        ('PEPOCH',r'Refrence epoch for $\alpha$, $\delta$, and $\nu$ (MJD)'),
        ('PMRA',r'Proper motion in $\alpha$, $\nu_{\alpha}=\dot{\alpha}\cos \delta$ (mas~yr$-1$)'),
        ('PMDEC',r'Proper motion in $\delta$, $\nu_{\delta}=\dot{\delta}$ (mas~yr$-1$)'),
        ('PX',r'Parallax, $\pi$ (mas)'),
        ('DM',r'Dispersion Measure# (pc~cm$^{-3}$)'),
        ('SINI', r'$\sin i$, where $i$ is the magnetic inclination angle'),
        ('PB',r'Orbital Period, $P_{\rm b}$ (day)'),
        ('ECC',r'Eccentricity, $e$'),
        ('OM',r'Angle of periastron#, $\omega$ (deg)'),
        ('T0',r'Time of periastron passage, $T_0$ (MJD)'),
        ('PBDOT',r'Change rate of $P_{\rm b}$, $\dot{P}_{\rm b}$ ($10^{-12}$s~s$^{-1}$)'),
        ('A1',r'Projected semi-major axis, $x$ (lt-s)'),
        ('XDOT',r'Change rate of projected semi-major axis due to proper motion $\dot{x}$ (lt-s~s$^{-1}$)'),
        ('M2',r'Companion Mass, $M_c$ ($M_{\odot}$)'),
        ('FD1',r'Profile frequency dependency parameter, FD1 '),
        ('FD2',r'Profile frequency dependency parameter, FD2 '),
        ('FD3',r'Profile frequency dependency parameter, FD3 '),
        ('FD4',r'Profile frequency dependency parameter, FD4 ')
        ])

Fixed = OrderedDict([
    ('EPHEM', r'Solar System Ephemeris'),
    ('PEPOCH', r'Reference epoch (MJD)')
    ])

Derived = OrderedDict([
        ('KIN',r'Oribital inclination, $i$ (deg)'),
        ('KOM',r'Position angle of ascending node, $\Omega$ (deg)'),
        ('M1',r'Pulsar mass, $M_{\rm PSR}$ ($M_{\odot}$)'),
        ('B',r'Dipole magnetic field, $B$ (G)'),
        ('Age',r'Characteristic age, $\tau_c$ (yr)'),
        ])


Caption = 'Timing model parameters#.'
colnames = ['Parameter', 'Value']
parameter = []
value = []
data = []

parameter.append(r'\textit{Measured Parameters}')
value.append('')

for k in Observable:
    parameter.append(Observable[k])
    if type(m.__dict__[k]) in (tuple, list):
        val, err = m.__dict__[k]
        #if not all([type(v) is decimal.Decimal for v in m.__dict__[k]]):
        if k == 'RAJ':
            #print k, m.__dict__[k], [type(v) is decimal.Decimal for v in m.__dict__[k]]
            ra = RA(m.__dict__[k][0])
            raerr = m.__dict__[k][1]
            value.append('%s:%s:%s' % (ra.HH,ra.MM, parseerror(Decimal(ra.SS),Decimal(raerr))))
        elif k == 'DECJ':
            dec = Dec(m.__dict__[k][0])
            decerr = m.__dict__[k][1]
            value.append('%s:%s:%s' % (dec.dd,dec.mm, parseerror(Decimal(dec.ss),Decimal(decerr))))
        else:
            value.append(parseerror(*m.__dict__[k]))
    else:
        if k == 'DM':
            value.append(parseerror(*aveDM(m)))
        else:
            value.append((m.__dict__[k]))

parameter.append(r'\textit{Fixed Parameters}')
value.append('')

for k in Fixed:
    parameter.append(Fixed[k])
    try:
        value.append(parseerror(*m.__dict__[k]))
    except:
        value.append(SF(globals()[k](m)))
data = [parameter, value]

parameter.append(r'\textit{Derived Parameters}')
value.append('')

for k in Derived:
    parameter.append(Derived[k])
    try:
        value.append(parseerror(*m.__dict__[k]))
    except:
        value.append(SF(globals()[k](m)))
data = [parameter, value]


table = deluxetable(Caption=Caption, colsetting='lc', colnames = colnames, data=data, label="tab:par", comments=['Numbers in parentheses indicate the uncertainties on the last disgit(s).  Uncertainties on parameters are estimated using MCMC simulation. ', 'The averaged DM value; See Section \label{sec:dmx} and Figure \label{fig:dmx} for the more discussion.', 'See Figure 2 of \citealt{sns+05} for definition.'])

print table


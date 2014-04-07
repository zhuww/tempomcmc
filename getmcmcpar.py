from tempo import *
import cPickle as pickle


pf = model('bestmcmc.par')

dit = pickle.load(open('bestpar.p', 'rb'))
best = dit['BEST']
plist = dit['parameters'] 
MChain = pickle.load(open('MChain.p','rb'))
MarkovChain = MChain['Chain']
print len(MarkovChain)
pi = 3.141592653589793
G = 6.67384e-11
Msun = 1.98892e30
c = 2.99792458e8
twopi = 6.283185307179586
fac = 1.536e-16 
Tsun = 4.925490947

pf.ECC = pf.E
pf.parameters['ECC'] = '1'
pf.manifest.append('ECC')

for k in plist:
    if k in ['RAJ' , 'DECJ' ]:
        idx = plist.index(k)
        firstval = ':'.join(MarkovChain[0][idx].split(':')[:-1])
        val = np.array([float(p[idx].split(':')[-1]) for p in MarkovChain])
        m = sum(val)/len(val) 
        var = [float(x) for x in val - m]
        e = np.std(var)
        ss = str(m)
        ds, fs = ss.split('.')
        pf.__dict__[k] = (firstval+':'+ds.zfill(2)+'.'+fs, Decimal(str(e)))

    #elif k in ['F0', 'F1', 'T0']:
    else:
        idx = plist.index(k)
        val = np.array([(p[idx]) for p in MarkovChain])
        m = sum(val)/len(val) 
        var = [float(x) for x in val - m]
        e = np.std(var)
        pf.__dict__[k] = (m, Decimal(str(e)))
    #else:
        #idx = plist.index(k)
        #val = np.array([float(p[idx]) for p in MarkovChain])
        #pf.__dict__[k] = (Decimal(str(val.mean())), Decimal(str(np.std(val))))

im2 = plist.index('M2')
ipb = plist.index('PB')
isini = plist.index('SINI')
ia = plist.index('A1')
ichisq = plist.index('chisq')
M2 = np.array([float(p[im2]) for p in MarkovChain])
Pb = np.array([float(p[ipb])*secperday for p in MarkovChain])
SINI = np.array([float(p[isini]) for p in MarkovChain])
a = np.array([float(p[ia])*c for p in MarkovChain])
#M1 = (Pb/2/np.pi*np.sqrt(G*(M2*SINI)**3/a**3)-M2)/Msun
M1 = (Pb/2/np.pi*np.sqrt(Tsun*(M2*SINI)**3/a**3)-M2)
#M2 = M2/Msun
KIN = np.arcsin(SINI)*180./np.pi
chisq = [p[ichisq] for p in MarkovChain]
bestidx = chisq.index(min(chisq))
k = 'M1'
pf.__dict__[k] = (Decimal(str(M1.mean())), Decimal(str(np.std(M1))))
pf.parameters[k] = '1'
pf.manifest.append(k)
k = 'KIN'
pf.__dict__[k] = (Decimal(str(KIN.mean())), Decimal(str(np.std(KIN))))
pf.parameters[k] = '1'
pf.manifest.append(k)

#for p in set(plist) & set(pf.parameters):
    #idp = plist.index(p)
    #pf.__dict__[p][0] = best[idp]

if 'PAASCNODE' in plist:
    pf.KOM = pf.PAASCNODE
    pf.parameters['KOM'] = '1'
    pf.manifest.append('KOM')
    p = 'PAASCNODE'
    idp = plist.index(p)
    pf.__dict__[p] = best[idp]

pf.write('mcmcresult.par')

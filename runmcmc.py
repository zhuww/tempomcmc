"""
An code for running MCMC simulation to determine the confidence range of tempo parameters. Code in use:
    runmcmc.py : the main driving program
    plotmc.py : the plotting program
    ProgressBar.py : for plotting the progress bar
"""
from tempo import tempofit, tempo2fit, touchparfile, uniquename, PARfile #, model, TOAfile
from math import *
from decimal import *
import os
from copy import *
from numpy.random import normal , uniform ,seed
import numpy as np
import time

from multiprocessing import Pool, Process, Manager

def randomnew(pf, stepsize): #special for 1713
    twopi = 6.283185307179586
    fac = 1.536e-16 * 1.e12
    x = float(str(pf.A1[0]))
    sini = float(str(pf.SINI[0]))
    cosi = -1. * np.sqrt(1 - sini**2)
    Omega = float(str(pf.PAASCNODE))
    m2 = float(str(pf.M2[0])) + normal(0,0.05*stepsize)
    cosi = cosi + uniform(-0.5, 0.5, 1)[0]*0.001*stepsize #normal(0, 0.001*stepsize)
    Omega = Omega + normal(0, 4.0*stepsize)
    mu = np.sqrt(float(str(pf.PMRA[0]**2+pf.PMDEC[0]**2)))
    #print 'mu:', mu
    #sini = sqrt(1 - cosi**2)
    thetamu = 180. + np.arctan(float(str(pf.PMRA[0]/pf.PMDEC[0])))/np.pi*180
    xdot = -1.* fac * x * mu * (cosi/sini) * sin((thetamu-Omega)*twopi/360.)
    sini = np.sqrt(1 - cosi**2)
    pf.SINI[0] = Decimal(str(sini))
    pf.XDOT[0] = Decimal(str(xdot))
    pf.PAASCNODE = Decimal(str(Omega))
    pf.M2[0] = Decimal(str(m2))
    #print np.arcsin(sini)*180/np.pi, Omega, thetamu, xdot
    return pf

def probcal(pf):
    global smallestchisq
    pf.write()
    #m2 = float(str(pf.M2[0]))
    #Omega = float(str(pf.PAASCNODE))
    #sini = float(str(pf.SINI[0]))
    #if m2 <= 0 or Omega > 360 or Omega < -360 or sini > 1.:
        #return 0
    chisq, dof = tempofit(parfile, toafile = toafile, pulsefile = pulsefile)
    pf.chisq = chisq
    print dof, chisq
    if chisq < smallestchisq: smallestchisq = chisq
    try:
        return exp((smallestchisq - chisq)/2.) #Ingrid/Paul?
    except OverflowError:
        print chisq, smallestchisq
        print pf.parfile
        raise OverflowError 

#print probcal(90, 0.25, 0.28)


from itertools import count
import cPickle as pickle
import os,sys
from tempfile import mkdtemp

#print probcal(iOmega, icosi, im2)

manager = Manager()

class MChain(object):
    def __enter__(self):
        #try:
            #Chain = pickle.load(open('MChain.p', 'r'))['Chain']
        #except:
            #Chain = []
        self.Chain = []
        self.cwd = os.getcwd()
        return self
    def __exit__(self, exc_type, exc_value, exc_tb):
        os.chdir(self.cwd)
        #try:
            #f = open(self.cwd+'/MChain.p', 'rb')
            #MarkovChain = pickle.load(f)['Chain']
            #f.close()
            #f1 = open(self.cwd+'/MChain.p'+str(os.getpid()), 'rb')
            #ThisChain = pickle.load(f1)['Chain']
            #f1.close()
            #MarkovChain.extend(ThisChain)
            #MarkovChain.extend(self.Chain)
        #except IOError:
            #MarkovChain = self.Chain
        #except EOFError:
            #print 'encouter EOFerror at exit'
            #time.sleep(10)
            #self.__exit__(exc_type, exc_value, exc_tb)
            #return True
        #if len(MarkovChain)>2: 
            #if len(MarkovChain[-1]) < len(MarkovChain[-2]):
                #MarkovChain = MarkovChain[:-1]
        #dit = {'Chain':MarkovChain}
        #f = open(self.cwd+'/MChain.p', 'wb')
        #pickle.dump(dit, f, protocol=2)
        #f.flush()
        #f.close()
        #del dit
        #os.remove(self.cwd+'/MChain.p'+str(os.getpid()))
        #print len(MarkovChain), 'points saved to MChain.p'
        #self.save()

        try:
            os.remove(self.cwd+'/MChain.p'+str(os.getpid()))
        except:pass
        if exc_type is KeyboardInterrupt:
            print '\nManually Stopped\n'
            return True
        else:
            return exc_type is None
        print '\nFinish running\n' 

    def save(self):
        try:
            #f = open(self.cwd+'/MChain.p', 'rb')
            #MarkovChain = pickle.load(f)['Chain']
            #f.close()
            MarkovChain.extend(self.Chain)
            self.Chain = []
        except IOError:
            MarkovChain = self.Chain
        except EOFError:
            print 'encounter EOFError here', os.getpid(), cwd
            #time.sleep(10)
            self.save()
            return
        if len(MarkovChain)>2: 
            if len(MarkovChain[-1]) < len(MarkovChain[-2]):
                MarkovChain = MarkovChain[:-1]
        dit = {'Chain':np.array(MarkovChain)}
        f = open(self.cwd+'/MChain.p'+str(os.getpid()), 'wb', 0)
        pickle.dump(dit, f, protocol=2)
        f.flush()
        f.close()
        del dit
        #print 'saved per 100 points',  len(MarkovChain)

def motifile(file, cwd, tmpdir):
    os.system('cp %s/%s %s/%s' % (cwd, file, tmpdir, file))
    text = ''
    f = open(file, 'rw')
    for l in f.readlines():
        if not l.find('INCLUDE') == -1:
            a = l.split()
            if a[0] == 'C' or a[0] =='#':
                continue
            if not open(cwd+'/'+a[1],'r').read().find('INCLUDE') == -1: 
                motifile(a[1], '..', '.')
                l = a[0] +' '+a[1]
            else:
                l = a[0] + ' '+cwd+'/'+a[1]
            if not l[-1] == '\n':
                l += '\n'
            text += l
        else:
            if not l[-1] == '\n':
                l += '\n'
            text += l
    f.close()
    f = open(file, 'w')
    f.write(text)
    f.close() #motify the tim file to make sure INCLUDE follow the right files.
    
from ProgressBar import progressBar
def mcmc(Chain, runtime, MarkovChain, mixingtime=1000, stepsize=1, seed=0 ):
    pb = progressBar(maxValue = runtime + mixingtime)
    cwd=os.getcwd()
    tmpdir = cwd+'/.'+uniquename()
    if not tmpdir == None:
        if os.path.exists(tmpdir):
            os.chdir(tmpdir)
        else:
            os.mkdir(tmpdir)
            os.chdir(tmpdir)
    os.system('cp %s/%s %s/%s' % (cwd, parfile, tmpdir, parfile))
    os.system('cp %s/%s %s/%s' % (cwd, pulsefile, tmpdir, pulsefile))
    motifile(toafile, cwd, tmpdir)
    touchparfile(parfile, NITS=1)
    pf = PARfile(parfile)
    #chisq, dof = tempofit(parfile, toafile = toafile, pulsefile = pulsefile)
    pf.matrix(toafile)
    pf.freezeall()
    #pf.thawall('JUMP_')
    #pf.write()
    #plist = [x for x in pf.manifest if x in pf.parameters.keys() ]
    #dit = {'BEST':[pf.__dict__[p][0] for p in plist] + [ chisq], 'parfile':pf.parfile, 'parameters':plist + [ 'chisq']}
    #pickle.dump(dit, open('%s/bestpar.p' % cwd, 'w', 0), protocol=2)
    p0 = probcal(pf)
    p = p0
    ThisChain = []
    c = 0
    while c <= mixingtime + runtime:
        c+=1
        npf = pf.randomnew(stepsize=stepsize)
        #randomnew(npf, stepsize) #only use this for 1713
        p1 = probcal(npf)
        if c % 30 == 0:pb(c)
        if c > mixingtime and c % (100+(seed%100)) == 0:
            MarkovChain.extend(Chain.Chain)
            ThisChain.extend(Chain.Chain)
            Chain.Chain = [] #empty the list
            #try:
            TC = np.array(ThisChain)
            dit = {'Chain':TC}
            pid = str(os.getpid())
            try: 
                os.remove(cwd+'/MChain.p'+pid)
            except:pass
            f = open(cwd+'/MChain.p'+pid, 'wb', 0)
            pickle.dump(dit, f, protocol=2)
            f.flush()
            f.close()
            del dit
            del TC
        if p1 > p0:
            if c > mixingtime:
                Chain.Chain.append([npf.__dict__[p][0] for p in plist] + [ npf.chisq])
            pf = npf
            p0 = p1
            if p1 > p:
                p = p1
                #if 'PAASCNODE' in plist:
                    #dict['BEST'] = [pf.__dict__[p][0] for p in plist[:-1]] + [pf.__dict__[p] for p in plist[-1:]] + [npf.chisq]
                #else:
                dit['BEST'] = [npf.__dict__[p][0] for p in plist] + [ npf.chisq]
                pickle.dump(dit, open('%s/bestpar.p' % cwd, 'wb', 0), protocol=2)
        else:
            t = uniform(0,1,1)[0]
            if t < p1/p0:
                if c > mixingtime:
                    Chain.Chain.append([npf.__dict__[p][0] for p in plist] + [ npf.chisq])
                #print npf.M2[0], npf.chisq
                pf = npf
                p0 = p1
            else:
                if c > mixingtime:
                    Chain.Chain.append([pf.__dict__[p][0] for p in plist] + [ npf.chisq])
    MarkovChain.extend(Chain.Chain)
    os.chdir(cwd)

        
from optparse import OptionParser
if __name__ == '__main__':
    #main()
    usage = "usage: %prog [options] arg"
    parser = OptionParser()
    parser.add_option("-f", '--parfile', dest="parfile", help="par file")
    parser.add_option("-t", '--timfile', dest="toafile", help="toa file")
    parser.add_option("-n", '--pulsefile', dest="pulsefile", help="pulse number file", default=None)
    parser.add_option("-i", '--iter', type='int', nargs=1, dest='steps', help="number of steps")
    parser.add_option("-m", '--mixing', type='int', nargs=1, dest='mixing', help="number of mixing steps", default=1000)
    parser.add_option("-p", '--parallel', type='int', nargs=1, dest='paral', help="number of parallel processes")
    parser.add_option("-s", '--seed', type='int', nargs=1, dest='seed', default=int(os.getpid()), help="random number seed")
    parser.add_option("-z", '--stepsize', type='float', nargs=1, dest='stepsize', default=1., help="step size")
    (options, args) = parser.parse_args(args=sys.argv[1:])
    print options

    parfile = options.parfile
    toafile = options.toafile
    pulsefile = options.pulsefile
    steps = options.steps
    mixing = options.mixing
    rseed = options.seed
    stepsize = options.stepsize
    px = options.paral
    pf =PARfile(parfile)
    #pf = model(parfile)
    pf.freezeall()
    pf.thawall('JUMP_')
    pf.write('mcmc.par')
    touchparfile('mcmc.par', NITS=1)
    chisq, dof = tempofit('mcmc.par', toafile = toafile, pulsefile = pulsefile)
    #pf.tempofit(TOAfile(toafile), pulsefile = pulsefile)
    smallestchisq = chisq
    #print 'smallestchisq', smallestchisq, dof

    cwd=os.getcwd()
    plist = [x for x in pf.manifest if x in pf.parameters.keys() if not x.startswith('DMX') and not x.startswith('JUMP')]

    dit = {'BEST':[pf.__dict__[p][0] for p in plist] + [chisq], 'parfile':pf.parfile, 'parameters':plist + ['chisq']}
    pickle.dump(dit, open('%s/bestpar.p' % cwd, 'w', 0), protocol=2)
    



    def run(argv):
        s, MarkovChain = argv
        seed(s) # assigning different initial seed for the random number generator in different threads.
        #steps = 50000
        with MChain() as Chain:
            #print type(Chain)
            #print steps, tmpdir
            mcmc(Chain, steps, MarkovChain, mixingtime=mixing, stepsize=stepsize, seed=s)
        return Chain

    MarkovChain = manager.list()
    cwd = os.getcwd()

    if px == None:
        run([rseed, MarkovChain])
    else:
        px = options.paral
        p = Pool(px)
        p.map(run, [(rseed+100*s, MarkovChain) for s in range(px)])

    #print 'generated %s pints' % len(MarkovChain)
    try:
        f = open('./MChain.p', 'rb')
        oldChain = list(pickle.load(f)['Chain'])
        f.close()
        #oldChain.extend(MarkovChain)
        MarkovChain = oldChain + list(MarkovChain)
    except IOError:pass
    dit = {'Chain':np.array(MarkovChain)}
    f = open(cwd+'/MChain.p', 'wb')
    pickle.dump(dit, f, protocol=2)
    #f.flush()
    f.close()





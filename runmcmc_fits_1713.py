from math import *
from tempo import tempofit, tempo2fit, touchparfile, uniquename, PARfile #, model, TOAfile
from numpy.random import normal , uniform ,seed
import numpy as np
import time
from multiprocessing import Pool, Process, Manager, JoinableQueue, Value
import os,sys
import pyfits
import fitsio
from itertools import count
import cPickle as pickle
import os,sys
from tempfile import mkdtemp
from decimal import *
from copy import deepcopy


def randomnew(pf, stepsize): #special for 1713
    twopi = 6.283185307179586
    fac = 1.536e-16 * 1.e12
    x = float(str(pf.A1[0]))
    sini = float(str(pf.SINI[0]))
    cosi = -1. * np.sqrt(1 - sini**2)
    Omega = float(str(pf.PAASCNODE))
    m2 = float(str(pf.M2[0])) + normal(0,0.01*stepsize)
    cosi = cosi + uniform(-0.5, 0.5, 1)[0]*0.001*stepsize #normal(0, 0.001*stepsize)
    Omega = Omega + normal(0, 0.1*stepsize)
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
    #if chisq < smallestchisq: smallestchisq = chisq
    try:
        #return exp((smallestchisq - chisq)/2.) #Ingrid/Paul?
        return (smallestchisq - chisq)/2. #Ingrid/Paul?
    except OverflowError:
        print chisq, smallestchisq
        print pf.parfile
        raise OverflowError 

#print probcal(90, 0.25, 0.28)

def initfitsfile(pf):
    plist = [x for x in pf.manifest if x in pf.parameters.keys() if not x.startswith('DMX') and not x.startswith('JUMP')] + ['PAASCNODE']
    cols = []
    for par in plist:
        if par in ['RAJ', 'DECJ', 'RA', 'DEC']:
            col = pyfits.Column(name=par, format='D', array=[])
        elif par in pf.LongParameters:
            col = pyfits.Column(name=par, format='D', array=[])
        else:
            #col = pyfits.Column(name=par, format='D', array=[float(pf.__dict__[par][0])])
            col = pyfits.Column(name=par, format='D', array=[])
        cols.append(col)
    #cols.append(pyfits.Column(name='chisq', format='D', array=[chisq]))
    cols.append(pyfits.Column(name='chisq', format='D', array=[]))
    newtbl = pyfits.new_table(cols)
    newhdr = newtbl.header
    for par in [p for p in plist if p in ['RAJ', 'DECJ', 'RA', 'DEC'] or p in pf.LongParameters]:
        newhdr.set(par , str(pf.__dict__[par][0]))
    newhdr.set('SPCPAR', '|'.join(pf.LongParameters))
    PSRname = pf.PSR
    try:
        newtbl.writeto(PSRname+'.mcmc')
    except IOError:
        print 'file %s already exist, would you like to append to the existing file?(y/n)' % (PSRname+'.mcmc')
        if raw_input().lower() == 'y':
            pass
        else:
            sys.exit(0)



class MChain(object):
    def __enter__(self):
        self.Chain = []
        self.cwd = os.getcwd()
        return self
    def __exit__(self, exc_type, exc_value, exc_tb):
        os.chdir(self.cwd)

        #try:
            #os.remove(self.cwd+'/MChain.p'+str(os.getpid()))
        #except:pass
        if exc_type is KeyboardInterrupt:
            print '\nManually Stopped\n'
            return True
        else:
            return exc_type is None
        print '\nFinish running\n' 

    def save(self):
        try:
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

from ProgressBar import progressBar
def mcmc(Chain, runtime, MarkovChain, Global, mixingtime=1000, stepsize=1, seed=0 ):
    global smallestchisq
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
    pf0 = PARfile(parfile)
    pforig = PARfile(cwd+'/'+parfile)
    pf0.LongParameters = []
    pf = pf0
    for par in [p for p in pf.parameters if not p in ['RAJ', 'DECJ']]:
        val,err = pf.__dict__[par]
        if err > 0 and not pf.parameters[par] == 0:
            NumSigDig = np.log10(np.abs(val/err)) 
            if NumSigDig > 10:
                pf.LongParameters.append(par)
    plist = [x for x in pf.manifest if x in pf.parameters.keys() if not x.startswith('DMX') and not x.startswith('JUMP')] + ['PAASCNODE'] 
    dtypes = np.dtype([(p, '>f8') for p in (plist+['chisq'])])
    pf0.matrix(toafile)
    pf0.freezeall()
    p0 = probcal(pf0)
    Global.pmax = p0
    #print 'inital pmax', Global.pmax, p0+smallestchisq
    ThisChain = []
    c = 0
    randomlist = uniform(0,1,size=runtime+mixingtime)
    def savepar(npf, pf, plist):
        dataarray = []
        for p in plist:
            if p in ['RAJ', 'DECJ']:
                dataarray.append(float(npf.__dict__[p][0].split(':')[-1]) - float(pf.__dict__[p][0].split(':')[-1]))
            elif p in pf.LongParameters:
                dataarray.append(float(str(npf.__dict__[p][0] - pf.__dict__[p][0])))
            elif p == 'PAASCNODE':
                dataarray.append(float(npf.__dict__[p]))
            else:
                dataarray.append(float(npf.__dict__[p][0]))
        dataarray.append(npf.chisq)
        return tuple(dataarray)

    while c <= mixingtime + runtime - 1:
        c+=1
        npf = pf.randomnew(stepsize=stepsize)
        randomnew(npf, stepsize) #only use this for 1713
        p1 = probcal(npf)
        if c % 30 == 0:pb(c)
        t = randomlist[c-1]
        if t < exp(p1-p0):
            pf = npf
            p0 = p1
            if p1 > Global.pmax:
                Global.pmax = p1
                npf0 = deepcopy(pf)
                npf0.parameters = pforig.parameters
                npf0.write(cwd+'/'+PSRname+'.mcmcbest.par')
                print '\nnew best parfile saved to:', cwd+'/'+PSRname+'.mcmcbest.par'
                print 'pmax:', Global.pmax, 'new chisq:', pf.chisq, 'old chisq:', smallestchisq
        if c > mixingtime:
            if t < exp(p1-p0):
                Chain.Chain.append(savepar(npf, pf0, plist))
            else:
                Chain.Chain.append(savepar(pf, pf0, plist))
        if c % (100+(seed%100)) == 0:
            data = np.array(Chain.Chain, dtype=dtypes)
            MarkovChain.put(data)
            Chain.Chain = []
    data = np.array(Chain.Chain, dtype=dtypes)
    MarkovChain.put(data)

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

from optparse import OptionParser
if __name__ == '__main__':
    usage = "usage: %prog [options] arg"
    parser = OptionParser()
    parser.add_option("-f", '--parfile', dest="parfile", help="par file")
    parser.add_option("-t", '--timfile', dest="toafile", help="toa file")
    parser.add_option("-n", '--pulsefile', dest="pulsefile", help="pulse number file", default=None)
    parser.add_option("-i", '--iter', type='int', nargs=1, dest='steps', help="number of steps")
    parser.add_option("-m", '--mixing', type='int', nargs=1, dest='mixing', help="number of mixing steps", default=1000)
    parser.add_option("-p", '--parallel', type='int', nargs=1, dest='paral', help="number of parallel processes", default=1)
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
    pf = PARfile(parfile)
    pf.LongParameters = []
    for par in [p for p in pf.parameters if not p in ['RAJ', 'DECJ']]:
        val,err = pf.__dict__[par]
        if err > 0 and not pf.parameters[par] == 0:
            NumSigDig = np.log10(np.abs(val/err)) 
            if NumSigDig > 10:
                pf.LongParameters.append(par)
    pf.freezeall()
    pf.thawall('JUMP_')
    pf.write('mcmc.par')
    touchparfile('mcmc.par', NITS=1)
    chisq, dof = tempofit('mcmc.par', toafile = toafile, pulsefile = pulsefile)
    #pf.tempofit(TOAfile(toafile), pulsefile = pulsefile)
    smallestchisq = chisq
    initfitsfile(pf)
    PSRname = pf.PSR

    def save_chain(chain):
        ff = fitsio.FITS(PSRname+'.mcmc', 'rw')
        ff[-1].append(chain)
        ff[-1].write_checksum()
        ff.close()

    def collector(queue, savecount=10):
        chain = queue.get()
        counter = 1
        while True:
            newchain = queue.get()
            chain = np.append(chain, newchain) 
            counter += 1
            if counter >= savecount:
                #print 'save:', chain.shape
                save_chain(chain)
                chain = queue.get()
                counter = 1
            #queue.task_done()
        save_chain(chain)
        queue.task_done()

    class worker(Process):
        def __init__(self, MarkovChain, Global,  steps, mixing, stepsize, seed ):
            Process.__init__(self)
            self.queue = MarkovChain
            self.steps = steps
            self.mixing = mixing
            self.stepsize = stepsize
            self.seed = seed
            self.Global=Global 
        def run(self):
            #s, MarkovChain = argv
            np.random.seed(self.seed) # assigning different initial seed for the random number generator in different threads.
            with MChain() as Chain:
                mcmc(Chain, self.steps, self.queue, self.Global, mixingtime=self.mixing, stepsize=self.stepsize, seed=self.seed)
            self.queue.task_done()
            return 

    #MarkovChain = manager.list()
    MarkovChain = JoinableQueue()
    manager = Manager()
    Global  = manager.Namespace()
    Global.pmax = 0.
    cwd = os.getcwd()
    if px == None:
        ChainSaver = Process(target=collector, args=((MarkovChain,)))
    else:
        ChainSaver = Process(target=collector, args=((MarkovChain,px)))
    ChainSaver.daemon = True
    ChainSaver.start()

    rseed = int(os.getpid())
    if px == None:
        #run([rseed, MarkovChain])
        w = worker(MarkovChain, Global, steps, mixing, stepsize, rseed)
        w.start()
        MarkovChain.join()
        w.join()
    else:
        px = options.paral
        #p = Pool(px)
        workforce = []
        for s in range(px):
            seed = rseed+100*s
            workmule = worker(MarkovChain, Global, steps, mixing, stepsize, seed)
            workforce.append(workmule)
        for w in workforce:
            w.start()
        #p.map(run, [(rseed+100*s, MarkovChain) for s in range(px)])
        MarkovChain.join()
        for w in workforce:
            w.join()


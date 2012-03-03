'''A module that handles TOA files. 
'''
#import numpy, os
#from numpy import mean, std, var, array, double, angle
from math import *
from fileio import *
#from SaveLoadable import MetaSaveLoader
#from decimal import *
#from numpy import float64 as __Decimal
from decimal import Decimal as __Decimal
def Decimal(value):
    """Convert string segment from TOA line to a float number.
    return 0 when the input is an empty string.
    """
    if isinstance(value, (int,float)):
        return __Decimal(repr(value))
    elif isinstance(value, (basestring)):
        if value.strip() == '':
            return  __Decimal('0')
        else:
            return __Decimal(value)
    else:
        return __Decimal(value)


Observatory_list = {
        '1':   'GBT',
        '2':   'QUABBIN',
        '3':   'ARECIBO',
        '4':   'Hobart, Tasm',
        '5':   'PRINCETON',
        '6':   'VLA',
        '7':   'PARKES',
        '8':   'JODRELL BANK',
        '9':   'GB 300FT',
        'A':   'GB 140FT',
        'B':   'GB 85-3',
        'C':   'VLA SITE',
        'D':   'NORTHERN CRO',
        'E':   'MOST',
        'F':   'Nancay',
        'G':   'Effelsberg',
        'H':   'JODRELL BANK',
        'I':   'JODRELL BANK',
        'J':   'mkiii',
        'K':   'tabley',
        'L':   'darnhall',
        'M':   'knockin',
        'N':   'defford',
        'P':   'JODRELL BANK',
        'Q':   'JODRELL BANK',
        'R':   'GMRT',
        'S':   'WSRT',
        'T':   'VLT',
        'U':   'VLT',
} 



class TOA(object):
    """A class for reading/formating a single TOA line """
    #__metaclass__ = MetaSaveLoader
    #reserved_words = ['MODE', 'JUMP']
    def __init__(self, line, file=None, **flags): #this way the keys in flags are passed as immutable values (parameters) to TOA, no need to do the copy/deepcopy.
        self.flags = flags
        if file:
            self.file = file
        else:
            self.file = ''
        if len(line) < 25:
            raise "Line too short."
        elif line[0] == 'C':
            raise "Commentted line."
        if line[0] == ' ':
            self.format = 'Parkes'
            self.info = []
            self.line = line.ljust(80, ' ')
            self._Observatory = self.line[79]
            self.Observatory = Observatory_list[str(self._Observatory)]
            self.info.append(self.line[1:25].strip())
            self.frequency = Decimal(self.line[25:34])
            self.TOA = Decimal(self.line[34:55])
            self.DMcorr = Decimal(self.line[55:63])
            self.TOAsigma = Decimal(self.line[63:71])
            self.info.append(self.line[71:80].strip())
        elif line[1] == ' ' and not line[0] == ' ':
            self.format = 'Princeton'
            self.info = []
            self.line = line.ljust(78, ' ')
            self._Observatory = int(self.line[0])
            self.Observatory = Observatory_list[str(self._Observatory)]
            self.info.append(self.line[2:15].strip())
            self.frequency = Decimal(self.line[15:24])
            self.TOA = Decimal(self.line[24:44])
            self.TOAsigma = Decimal(self.line[44:53])
            self.info.append(self.line[53:68].strip())
            self.DMcorr = Decimal(self.line[68:77])
        elif not line[0] ==' ' and not line[1] == ' ' and line[14] == '.':
            self.format = 'ITOA'
            self.info = []
            self.line = line.ljust(59, ' ')
            self._Observatory = self.line[57:59]
            self.Observatory = Observatory_list[str(self._Observatory.strip())[0]]
            self.info.append(self.line[:9].strip())
            self.TOA = Decimal(self.line[9:28])
            self.TOAsigma = Decimal(self.line[28:34])
            self.frequency = Decimal(self.line[34:45])
            self.DMcorr = Decimal(self.line[45:55])
            self.info.append(self.line[55:57].strip())
        else:
            self.format = 'Tempo2'
            self.info = []
            self.line = line
            line = self.line.split()
            self.info.append(line[0])
            self.frequency = Decimal(line[1])
            self.TOA = Decimal(line[2])
            self.TOAsigma = Decimal(line[3])
            self._Observatory = line[4]
            self.Observatory = Observatory_list[str(self._Observatory.strip())[0]]
            #self.info.extend(line[5:])
            tags = ' '.join(line[5:])
            import re
            p = re.compile('-(?P<par>\w+)\s+(?P<value>\S+)', re.VERBOSE)
            for m in p.finditer(tags):
                par = str(m.group('par'))
                value = str(m.group('value'))
                self.flags[par] = value


        #if self.flags.has_key('EQUAD'):
            #self.TOAsigma = sqrt(self.TOAsigma**2 + self.flags['EQUAD']**2)

    def __repr__(self):
        return '%s format TOA %s(%s us) at %s MHZ' % (self.format, self.TOA, self.TOAsigma, self.frequency)
    def __str__(self):
        if self.format == 'Princeton':
            #return '%d %s%9.3f%20.13f%9.3f%s%9.5f' % (self._Observatory, self.info[0][:13].ljust(13, ' '), self.frequency, self.TOA, self.TOAsigma,self.info[1][:15].rjust(15,' '), self.DMcorr)
            #return '%d %s%9.3f%s%9.3f%s%9.5f' % (self._Observatory, self.info[0][:13].ljust(13, ' '), self.frequency, str(self.TOA.quantize(Decimal(0.0000000000001))).rjust(20, ' '), self.TOAsigma,self.info[1][:15].rjust(15,' '), self.DMcorr)
            #return '%d %s%9.3f%s%9.3f%s%9.5f' % (self._Observatory, self.info[0][:13].ljust(13, ' '), self.frequency, str(self.TOA.quantize(Decimal(0.00000000000001))).rjust(20, ' '), self.TOAsigma,self.info[1][:15].rjust(15,' '), self.DMcorr)
            return '%d %s%8.2f %s%9.3f%s%9.5f' % (self._Observatory, self.info[0][:13].ljust(13, ' '), self.frequency, str(self.TOA.quantize(Decimal(0.00000000000001))).rjust(20, ' '), self.TOAsigma, self.info[1][:15].rjust(15,' '), self.DMcorr)
        if self.format == 'Parkes':
            return ' %s%8.2f %s%8.5f%8.3f%s%d' % (self.info[0].ljust(24, ' '), self.frequency,str(self.TOA.quantize(Decimal(0.0000000000001))).rjust(21, ' '), self.DMcorr, self.TOAsigma, self.info[1][:8].rjust(8, ' '), self._Observatory)
        if self.format == 'ITOA':
            try:
                info = ' '.join(self.info).strip().ljust(9, ' ')
            except TypeError:
                info = ' '.join(self.info).strip()[:9]

            if info[0] == ' ' or info[1] == ' ': info = 'IT'+info[2:]
            return '%s%s%6.3f%11.3f%10.5f  %s' % (info,str(self.TOA.quantize(Decimal(0.0000000000001))).rjust(19, ' '), self.TOAsigma,self.frequency, self.DMcorr, self._Observatory)
        if self.format == 'Tempo2':
            kwpars = ''
            for key in [k for k in sorted(self.flags.keys(), reverse=True) if not k == 'EQUAD' and not k == 'JUMPflag' and not k =='EMAX' and not k =='EFAC']:
                kwpars += ' -%s %s ' % (key, self.flags[key])
            return '%s %s %s %s %s %s' % (self.info[0], self.frequency, str(self.TOA.quantize(Decimal(0.0000000000001))).rjust(19, ' '), str(self.TOAsigma.quantize(Decimal(0.001))), self._Observatory, kwpars) 
    def tempo2fmt(self):
        file = self.file.replace(' ','_').replace('-', '_')
        fmtstr = '%s %s %s %s %s' % (file, self.frequency, str(self.TOA.quantize(Decimal(0.00000000000001))).rjust(20, ' '), self.TOAsigma, self._Observatory)
        kwpars = ''
        for key in [k for k in sorted(self.flags.keys(), reverse=True) if not k == 'EQUAD' and not k == 'JUMPflag' and not k =='EMAX' and not k == 'EFAC']:
            kwpars += ' -%s %s ' % (key, self.flags[key])
        return fmtstr+kwpars


class TOAcommand(object):
    def __init__(self, cmd, *args):
        self.cmd = cmd
        self.args = args
    def __str__(self):
        return '%s %s' % (self.cmd, ' '.join(self.args))
    def __repr__(self):
        return '%s %s' % (self.cmd, ' '.join(self.args))

class TOAfile(object):
    """[OutDated, use the TOAfile class in tempo module] A class for read/operate TOA files. """ 
    def __init__(self, file, kws={'JUMPflag':False}, F0=None): 
        self.list = []
        #self.PHASEJUMPS = {}
        self.EQUADvalues = {}
        self.start = Decimal('1000000')
        self.end = Decimal('0')
        self.F0 = F0
        with open(file,'rt') as f:
            text = f.read()
            lines = text.split('\n')
            toagroup = []
            SKIPflag = False
            for l in lines:
                s =l.strip().find('INCLUDE')
                if len(l) == 0:pass
                elif l[0] == 'C'  or l[0] == '#':
                    pass #commented
                elif not s == -1:
                    fn = l[s+7:].strip()
                    if not fn.find(' ') == -1:raise 'file name %s contains space?' % fn
                    tf = TOAfile(fn, kws=kws, F0=self.F0)
                    if tf.start < self.start: self.start = tf.start
                    if tf.end > self.end: self.end = tf.end
                    self.list.extend(tf.list)
                elif len(l) < 25:
                    if len(toagroup) > 0:self.list.append(toagroup)
                    command = TOAcommand(*l.split())
                    toagroup = []
                    if command.cmd == 'INFO':
                        INFO = command.args[0]
                        kws.update({'i': INFO})
                        self.list.append(command)#ignore
                    elif command.cmd == 'JUMP':
                        kws['JUMPflag'] = not kws['JUMPflag']
                        if kws['JUMPflag']:
                            if kws.has_key('jump'):
                                kws['jump'] += 1
                            else:
                                kws['jump'] = 1
                        #self.list.append(command)#ignore
                    elif command.cmd == 'PHASE':
                        PHASE = float(command.args[0])
                        #if F0 == None:
                            #print 'Need F0 to caculate phase offset.'
                            #raise Error
                        #Offset = Decimal(str(command.args[0]))/self.F0/86400
                        #if kws.has_key('p'):
                            #NewOffset = kws['p'] + Offset
                            #if NewOffset == 0:
                                #del kws['p']
                            #else:
                                #kws['p'] = NewOffset
                        #else:
                            #kws['p'] = Offset 
                        if kws.has_key('padd'):
                            NewPhase = kws['padd'] + PHASE
                            if NewPhase == 0.:
                                del kws['padd']
                            else:
                                kws['padd'] = NewPhase
                        else:
                            kws['padd'] = PHASE
                    elif command.cmd == 'EQUAD':
                        EQUAD = Decimal(command.args[0])
                        if not kws.has_key('f'):
                            key = 1
                            self.EQUADvalues[str(key)]=EQUAD
                        else:
                            key += 1
                            self.EQUADvalues[str(key)]=EQUAD
                        kws.update({'f':str(key)})
                    elif command.cmd == 'EMAX':
                        self.EMAX = Decimal(command.args[0])
                        kws['EMAX'] = self.EMAX
                        self.list.append(command)#ignore
                    elif command.cmd == 'EFAC':
                        self.EFAC = Decimal(command.args[0])
                        self.list.append(command)#ignore
                    elif command.cmd == 'SKIP':
                        SKIPflag = True
                    elif command.cmd == 'NOSKIP':
                        SKIPflag = False
                    else:
                        self.list.append(command)#ignore
                else:
                    if not SKIPflag:
                        if not file.find(' ') == -1:
                            file = ''
                        try:
                            t= TOA(l, file = file, **kws)
                        except:
                            print 'Probamatic line: ', l
                            print 'from file %s' % file
                        if self.__dict__.has_key('EMAX'): 
                            if self.__dict__.has_key('EFAC'):
                                if t.TOAsigma * self.EFAC > self.EMAX:
                                    continue 
                                else:toagroup.append(t)
                            else:
                                if t.TOAsigma > self.EMAX:
                                    continue
                                else:toagroup.append(t)
                        else:
                            toagroup.append(t)
                        if t.TOA < self.start:self.start = t.TOA
                        if t.TOA > self.end: self.end = t.TOA
                    else:
                        pass # skipped
            if len(toagroup) > 0:self.list.append(toagroup)
        def __str__(self):
            return '\n'.join([str(x) for x in self.list])
    def tempo2fmt(self):
        """convert TOA file to tempo2 format"""
        fmtstr = """
"""
        fmtstr += '\n'
        #toalist = [l for l in self.list if not isinstance(l,TOAcommand)]
        #toas = []
        #for l in toalist:
            #toas.extend(l)
        #toas.sort(key = lambda x:x.TOA)
        for toas in self.list:
            if isinstance(toas, (list,tuple)):
                fmtstr += '\n'.join([t.tempo2fmt() for t in toas])
                fmtstr += '\n'
            elif isinstance(toas, TOAcommand):
                fmtstr += '%s\n' % (toas)


        #for l in self.list:
            #if isinstance(l, TOAcommand):
                #fmtstr += str(l)
                #fmtstr += '\n'
            #else:
                #fmtstr += '\n'.join([t.tempo2fmt() for t in l])
                #fmtstr += '\n'
        #for key in self.PHASEJUMPS.keys():
            #fmtstr += 'PHASE %d %s 0\n' % (key, self.PHASEJUMPS[key])
        Extrastr = 'FORMAT 1\nMODE 1\n'
        for key in self.EQUADvalues.keys():
            Extrastr += 'T2EQUAD -f %s %s\n' % (key, self.EQUADvalues[key])
        return Extrastr + fmtstr
        #return '\n'.join([fmtstr] + self.list)


#from datatools.timing import ephemeris
def residual(file = 'resid2.tmp'):
    from numpy import fromfile
    data = fromfile(file='resid2.tmp', dtype=[
        ('N',           'i4'),
        ('toa',         'f8'),
        ('res_phase',   'f8'),
        ('res_sec',     'f8'),
        ('ophase',      'f8'),
        ('rf_bary',     'f8'),
        ('weight',      'f8'),
        ('err_us',      'f8'),
        ('prefit_sec',  'f8'),
        ('ddm',         'f8'),
        ('M',           'i4'),
        ])
    return data

# read jump offsets
def jumps(file = 'tempo.lis'):
    import re
    with open(file, 'r') as f: # open the file
        text = f.read()
        lines = text.split('\n')

    N = len(lines)

    p = re.compile('Iteration\s+(?P<interation>\d+)\s+of\s+(?P<total>\d+).*', re.VERBOSE) 
    e = re.compile('.*Weighted\s+RMS\s+residual.*', re.VERBOSE)

    # find out the useful part (the last interation)
    itr = 0
    istart = 0
    iend = N
    foundstart = False
    for i in range(N):
        l = lines[i]
        m = p.match(l) 
        if not m == None:
            nitr = m.group('interation') 
            #if nitr > itr:
                #itr = nitr
            istart = i
            foundstart = True
        if foundstart:
            if not e.match(l) == None:
                iend = i

    p = re.compile('PSR\s+\d{4}[-\+]\d{2,4}\s+Ephem.*', re.VERBOSE)  #find the parameter section
    foundstart = False
    for i in range(istart, iend):
        l = lines[i]
        if not p.match(l) == None:
            istart = i
            foundstart = True


    #write down the DMX parameters
    #DMX = {}
    #DMXmjdrange = {}
    #p = re.compile('(\s+(?P<par>(DM\s(Dot|Off))(\s*(?P<idx>\d+)))\s*)+', re.VERBOSE) # look for the DM parameter definition line. The values are 6 lines below the parameters.
    #e = re.compile('(?P<par>(DM\s(Dot|Off))(\s*(?P<idx>\d+)))', re.VERBOSE) 
    #for i in range(istart, iend):
        #l = lines[i]
        #if l.find('=') == -1 and l.find(':') == -1:
            #if not p.match(l) == None:
                #j = i
                #parlist = []
                #matchs = e.finditer(l)
                #for m in matchs:
                    #parlist.append(m.group('idx'))#, m.group('idx')
                #valuelist = lines[j+6].split()
                #rangelist = lines[j+2].split(' ')[1:]
                #for k in range(len(parlist)):
                    #DMX.update({parlist[k]:float(valuelist[k])})
                    #mjds, mjde = [float(x) for x in rangelist[k].split('-')]
                    #DMXmjdrange.update({parlist[k]: (mjds,mjde)})

    #maxidx = max([int(k) for k in DMX.keys()])
    #DMXarray = []
    #DMXrange = []
    #for i in range(maxidx):
        #DMXarray.append(DMX[str(i+1)])  # make sure no DMX parameter is missed
        #DMXrange.append(DMXmjdrange[str(i+1)])

    #search for and write down the offset parameters
    offset = {}
    mjdrange = {}
    p = re.compile('(\s+(?P<par>(Offset)(\s*(?P<idx>\d+)))\s*)+', re.VERBOSE) 
    e = re.compile('(?P<par>(Offset)(\s*(?P<idx>\d+)))', re.VERBOSE) 
    #idx = []
    for i in range(istart, iend):
        l = lines[i]
        if l.find('=') == -1 and l.find(':') == -1:
            if not p.match(l) == None:
                j = i
                parlist = []
                matchs = e.finditer(l)
                for m in matchs:
                    parlist.append(m.group('idx'))#, m.group('idx')
                valuelist = lines[j+6].split()
                rangelist = lines[j+2].split(' ')[1:]
                for k in range(len(parlist)):
                    offset.update({parlist[k]:float(valuelist[k])})
                    mjds,mjde = [float(x) for x in rangelist[k].split('-')]
                    mjdrange.update({parlist[k]: (mjds,mjde)})
    maxidx = max([int(k) for k in offset.keys()])
    offsetarray = []
    offsetrange = []
    for i in range(maxidx):
        offsetarray.append(offset[str(i+1)])  # make sure no offset parameter is missed
        offsetrange.append(mjdrange[str(i+1)])
    #inspect results
    from numpy import array
    #DMXarray = array(DMXarray)
    offarray = array(offsetarray)
    return offarray

import os,sys
from datetime import datetime

def uniquename():
    """Create an unique name for the working directory or the save files. The name start with when the python code was run %y-%m-%d_%Hh%Mm%Ss and then followed by '%s_%s_%s_%s' % (when,who,where,which), 'which' points to the process id.""" 
    who=os.environ['LOGNAME']
    where=os.uname()[1]
    which=os.getpid()
    when=datetime.now().strftime("%y-%m-%d_%Hh%Mm%Ss")
    what=__file__
    return "%s_%s_%s_%s" % (when,who,where,which)

def covtotempo2(timfile, parfile, newtimfile='', newparfile=''):
    """Run tempo and test if a tim file and a par works under tempo1, and then convert them to tempo2 format and add jumps to the new par file. """
    if not timfile or not parfile:
        raise 'must specify tim and par file'
    if os.access(newtimfile,os.R_OK):raise FileExist
    if os.access(newparfile,os.R_OK):raise FileExist
    try:
        cwd = os.getcwd()
        tempdir = '.' + uniquename()
        print tempdir
        os.mkdir(cwd+'/'+tempdir)
        os.chdir(cwd+'/'+tempdir)
        os.system('cp ../%s %s' % (timfile, timfile))
        os.system('cp ../%s %s' % (parfile, parfile))
        text = ''
        f = open(timfile, 'r')
        for l in f.readlines():
            if not l.find('INCLUDE') == -1:
                a = l.split()
                print a
                l = a[0] + ' ../'+a[1]
                if not l[-1] == '\n':
                    l += '\n'
                text += l
            else:
                if not l[-1] == '\n':
                    l += '\n'
                text += l
        f.close()
        f = open(timfile, 'w')
        f.write(text)
        f.close() #motify the tim file to enable those INCLUDE commands.

        os.system('tempo -f %s %s' % (parfile, timfile))
        resdata = residual()
        jps = jumps()
        #print len(resdata)
        print 'Number of Jumps: ', len(jps)
        #from datatools.timing import ephemeris
        parf = PARfile(parfile)
        #for i in range(len(jps)):
            #print 'Jump %d ' % i, '%f' % jps[i], parf.__dict__['JUMP_%d' % (i+1)]
        infogrp = []

        newparf = parf.psrname + '.par'
        f = open(newparf, 'a')
        i = 0
        for i in range(len(jps)):
            f.write('JUMP -jump %d' % (i+1) + ' %f 1\n' % jps[i])
        #for key in tf.PHASEJUMPS.keys():
            #f.write('PHASE %d %s 0\n' % (key, tf.PHASEJUMPS[key]))
        f.close()
        if newparfile:
            os.rename(newparf, newparfile)
            os.system('cp %s ../' % newparfile)
        F0 = PARfile(newparfile).F0[0]
        #print F0
        tf = TOAfile(timfile, F0=F0)
        for it in tf.list:
            if isinstance(it, TOAcommand):
                if it.cmd == 'INFO':
                    infogrp.append(it.args[0])
        if newtimfile:
            ntf = open(newtimfile, 'w')
            ntf.write(tf.tempo2fmt())
            ntf.close()
        os.system('cp %s ../' % newtimfile)
    finally:
        os.chdir(cwd)

from commands import getoutput
def tempofit(parfile, toafile=None, pulsefile=None):
    if toafile == None:
        toafile = ''
    if pulsefile == None:
        line = getoutput("tempo -w -c -f %s %s | grep 'Chisqr/nfree'" % (parfile, toafile))
    else:
        line = getoutput("tempo -w -c -f %s %s -ni %s | grep 'Chisqr/nfree'" % (parfile, toafile, pulsefile))
    #print line
    try:
        chisq = line.split('/')[1].split(':')[1]
    except:
        pass
        #print line
    try:
        dof = int(line.split('/')[2].strip(' ').split()[0])
    except:
        #print line
        dof =315
    try:
        chisq = float(chisq)
    except:
        chisq = 9999999.
    #print chisq, dof
    return chisq, dof

def tempo2fit(parfile, toafile=None):
    if toafile == None:
        toafile = ''
    line = getoutput("tempo2 -f %s %s -tempo1 > tmp2log; cat tmp2log | grep 'Chisqr/nfree'" % (parfile, toafile))
    chisqlist = []
    for l in line.split('\n'):
        l  = l.split()[6]
        chisq = l.split('/')[0]
        dof = int(l.split('/')[1])
        try:
            chisq = float(chisq)
        except:
            chisq = 99999.
        chisqlist.append(chisq)
    #print chisq, dof
    return min(chisqlist), dof

    
def touchparfile(parfile, **pars):
    parlist = [par for par in pars.keys()]
    with open(parfile, 'r') as parf:
        lines = parf.readlines()
        for i in range(len(lines)):
            try:
                key = lines[i].split()[0]
                if key in parlist:
                    if key == 'NITS':
                        lines[i] = '%s %s\n' % (key, pars[key])
                    else:
                        lines[i] = '%s %s 0\n' % (key, pars[key])
                    parlist.remove(key)
            except:pass
        for key in parlist:
            if key == 'NITS':
                lines.append('%s %s\n' % (key, pars[key]))
            else:
                lines.append('%s %s 0\n' % (key, pars[key]))
    with open(parfile, 'w') as parf:
        parf.write(''.join(lines))


from copy import deepcopy
from scipy.stats.kde import multivariate_normal as mvn
import numpy as np
from math import *
import decimal



params =[ '   f0','   f1','   f2','  Dec','   RA',
' pmdc',' pmra','    x','    e','   T0','   Pb','   Om',
' Omdt','gamma','   DM','   px',' Pbdt',' PPNg','    s',
'    M','   m2','  dth',' xdot',' edot','   x2','   e2',
'  T02','  Pb2','  Om2','   x3','   e3','  T03','  Pb3',
'  Om3',' PMRV',' XOmd',' Xpbd',' om2d','  x2d']

paramap = {
        '   f0':'F0',
        '   f1':'F1',
        '   f2':'F2',
        '  f03':'F3',
        '  f04':'F4',
        '  f05':'F5',
        '  f06':'F6',
        '  f07':'F7',
        '  f08':'F8',
        '  f09':'F9',
        '  Dec':'DECJ',
        '   RA':'RAJ',
        ' pmra':'PMRA',
        ' pmdc':'PMDEC',
        '    x':'A1',
        '    e':'E',
        '    s':'SINI',
        '   m2':'M2',
        '   T0':'T0',
        ' xdot':'XDOT',
        '   Pb':'PB',
        '   Om':'OM',
        ' Omdt':'OMDOT',
        ' Pbdt':'PBDOT',
        'gamma':'GAMMA',
        '   DM':'DM',
        '   px':'PX',
        ' O001':'JUMP_1',
        ' O002':'JUMP_2',
        ' O003':'JUMP_3',
        ' O004':'JUMP_4',
        ' O005':'JUMP_5',
        ' O006':'JUMP_6',
        ' O007':'JUMP_7',
        ' O008':'JUMP_8',
        ' O009':'JUMP_9',
        }
for i in range(100):
    key = 'DX' + str(i).rjust(3,'0')
    value = 'DMX_' + str(i).rjust(4,'0')
    paramap[key] = value
for i in range(100):
    key = 'D1' + str(i).rjust(3,'0')
    value = 'DMX1_' + str(i).rjust(4,'0')
    paramap[key] = value

class PARfile(object):
    #__metaclass__ = MetaSaveLoader
    def __init__(self, file):
        """A class object to contain the APIs for extracting ephemeris information from a .par file."""
        if os.access(file, os.R_OK):
            self.parfile = file
        else: 
            print file
            raise FileError(file) 
        def floatify(val):
            try: return Decimal(val)
            except: 
                import re
                p = re.compile('-*\d+\.*\d*D[-\+]\d+', re.VERBOSE)
                if p.match(val):
                    return Decimal(val.replace('D','E'))
                else:
                    return str(val)
        self.manifest = []
        self.parameters = {} # a place to hold the fit flag for each fitable parameter
        self.file = open(file, 'r')
        for lines in self.file.readlines():
            items = lines.split()
            if len(items) == 0:pass
            elif items[0] == '#':pass
            elif items[0][0] == '#':pass
            elif len(items) == 2 and not items[0] in ['SINI', 'M2', 'XDOT']: 
                self.__dict__[items[0]] = floatify(items[1])
                self.manifest.append(items[0])
            elif len(items) == 3 or items[0] in ['SINI', 'M2', 'XDOT']:
                value = floatify(items[1])
                if isinstance(value, decimal.Decimal):
                    error = value/10
                elif isinstance(value, basestring):
                    error = ''
                if items[0] in ['SINI', 'M2', 'XDOT']:
                    if len(items) == 2 and items[0] == 'SINI' and value == 'KIN':
                        self.__dict__[items[0]] = value
                        self.manifest.append(items[0])
                        continue
                    if len(items) < 4:
                        error = value/5
                        fitflag = '1' #set to fit for SINI, M2, and XDOT as default.
                    else:
                        fitflag = items[2]
                        error = floatify(items[3])
                else:
                    fitflag = items[2]
                    try:
                        if int(fitflag) == 1:
                            fitflag = '1'
                        else:
                            fitflag = '0'
                    except:
                        error = Decimal(fitflag)
                self.__dict__[items[0]] = [value, error]
                self.parameters.update({items[0]:fitflag})
                self.manifest.append(items[0])
            elif len(items) == 4:
                value = floatify(items[1])
                error = floatify(items[3])
                fitflag = items[2]
                self.parameters.update({items[0]:fitflag})
                self.__dict__[items[0]] = [value, error]
                self.manifest.append(items[0])
            elif len(items) == 5 and items[0] == 'JUMP':
                name = ' '.join(items[:3])
                self.manifest.append(name)
                value = floatify(items[3])
                fitflag = items[4]
                error = ''
                self.__dict__[name] = [value, error]
                self.parameters.update({name:fitflag})
            else:
                print items
                raise IndexError
            for key in self.parameters:
                try:
                    if int(self.parameters[key]) == 1:
                        self.parameters[key] = '1'
                    else:
                        self.parameters[key] = '0'
                except:
                    self.parameters[key] = '0'
        if self.__dict__.has_key('PSR'):
            self.psrname = self.PSR
        elif self.__dict__.has_key('PSRJ'): 
            self.psrname = 'J'+self.PSRJ
        elif self.__dict__.has_key('PSRB'): 
            self.psrname = 'B'+self.PSRB
        else:
            raise 'Cant Find PSR name'
        self.file.close()
    def write(self, fn=None):
        if fn == None:
            fn = self.parfile
        with open(fn, 'w') as f:
            text = ''
            for item in self.manifest:
                if item in self.parameters.keys():
                    text += '%s\t%s\t%s\t%s\n' % (item, self.__dict__[item][0], self.parameters[item] ,self.__dict__[item][1])
                else:
                    text += '%s\t%s\n' % (item, self.__dict__[item])
            f.write(text)
    def freezeall(self):
        for p in self.parameters:
            self.parameters[p] = '0'
    def thawall(self):
        for p in self.parameters:
            self.parameters[p] = '1'

    def matrix(self, timfile, TempoVersion=1):
        if self.BINARY == 'DDS':
            paramap['    s'] = 'SHAPMAX'
        if TempoVersion == 1:
            from numpy.core.records import fromfile
            os.system('tempo -j -f %s %s > tmplog' % (self.parfile, timfile))
            bpf = PARfile(self.PSR + '.par')
            par = []
            self.err = []
            cov = []
            self.parlist = []
            fd = open('matrix.tmp', 'rb')
            nn = fromfile(fd, formats='i4', shape=1)
            nn = int(nn[0][0])
            n = (nn - 29)/8

            for i in range(n):
                data = fromfile(fd, dtype=[
                    ('m',    '<i4'),
                    ('j',     '<i4'),
                    ('paramj','<a5'),
                    ('gcor',  '<f8'),
                    ('sig',   '<f8'), 
                    ], shape=1)
                j = int(data['j'][0])
                m = int(data['m'][0])
                a = fromfile(fd, formats='f8', shape=m, byteorder='<')
                cov.append([float(x[0]) for x in a])
                fd.seek(8, 1)
                paramj = str(data['paramj'][0])
                param = paramap[paramj]
                self.parlist.append(param)
                #try:
                self.err.append(float(str(bpf.__dict__[param][1])))
                #except:
                    #self.err.append(Decimal(0.05))
                #if not paramj in ['   RA', '  Dec']:
                    #par.append(float(str(pf.__dict__[paramap[paramj]][0])))
                #else:
                    #par.append(float(str(pf.__dict__[paramap[paramj]][0].split(':')[-1])))
                par.append(0)
            par = np.array(par)
            cov = np.matrix(cov)

            m,n = cov.shape
            if not m == n: raise exit
            #for i in range(m):
                #for j in range(n):
                    #cov[i,j] = cov[i, j] * err[i] * err[j]
            self.covariance = cov
            self.par = par
            Np = len(par)
            self.Nfac = int(sqrt(Np))+1
            #print self.Nfajc
        elif TempoVersion == 2:
            self.err = []
            os.system('tempo2 -f %s %s -tempo1 -newpar -output matrix > tempo2matrix.tmp' % (self.parfile, timfile))
            bpf = PARfile('new.par')
            with open('tempo2matrix.tmp' , 'r') as f:
                lines = f.readlines()
                i = 0
                for l in lines:
                    if not l.find('Correlation matrix') == -1:
                        break 
                    i += 1
                self.parlist = lines[i+2].split() 
                for par in self.parlist:
                    self.err.append(float(str(bpf.__dict__[par][1])))
                m = len(self.parlist)
                cov = []
                for j in range(m):
                    k = i+3+j
                    s = lines[k].split()[1:]
                    a = [float(x) for x in s]
                    a.extend([0]*(m-len(s)))
                    if not len(a) == m:
                        raise Error
                    cov.append(a)
                for i in range(m):
                    for j in range(i,m):
                        cov[i][j] = cov[j][i]
                self.covariance = np.matrix(cov)
                self.Nfac = (int(sqrt(m))+1)/3
                self.par = [0]*m
                #print self.covariance[:2,:2]
        else:
            raise TempoVersion
    def randomnew(self, stepsize=1.):
        new = deepcopy(self)
        err = mvn(self.par, self.covariance)
        for i in range(len(err)):
            p = self.parlist[i]
            if p == 'RAJ' or p == 'DECJ':
                dd,mm,ss = new.__dict__[p][0].split(':')
                #ss = str(Decimal(ss) + Decimal(str(err[i]*float(str(new.__dict__[p][1]))/self.Nfac)))
                #ss = str(Decimal(ss) + Decimal(repr(err[i]*self.err[i]/self.Nfac)))
                ss = str(Decimal(ss) + Decimal(repr(err[i]*self.err[i]/self.Nfac*stepsize)))
                new.__dict__[p][0] = ':'.join([dd,mm,ss])
            elif p == 'SINI':
                pass
            else:
                #new.__dict__[p][0] = new.__dict__[p][0] + Decimal(str(err[i]*float(str(new.__dict__[p][1]))/self.Nfac))
                new.__dict__[p][0] = new.__dict__[p][0] + Decimal(repr(err[i]*self.err[i]/self.Nfac*stepsize))

        return new


#get DMX parameters from a parfile
def getDMX(pf):

    #pf = PARfile('1713.all.par.t2')
    #pf2 = PARfile('newtempo2.zhuww.par')

    DMXpars = [p for p in pf.manifest if not p.find('DMX') == -1]
    DMX = {}
    DMXErr = {}
    DMX1= {}
    DMXR1 = {}
    DMXR2 = {}
    for pars in DMXpars:
        if pars.find('DMXR') == -1 and len(pars)>=4:
            if not pars[3] in 'EF':
                if not pars[3] == '1':
                    try:
                        DMX[int(pars[4:])] = pf.__dict__[pars][0] + pf.__dict__['DM']
                        DMXErr[int(pars[4:])] = pf.__dict__[pars][1] 
                    except:
                        print pars , pf.__dict__[pars]
                        raise stop
                else:
                    DMX1[int(pars[5:])] = pf.__dict__[pars][0]
        elif not pars.find('DMXR') == -1 and len(pars)>=4:
            if pars[4] == '1':
                DMXR1[int(pars[6:])] = pf.__dict__[pars]
            elif pars[4] == '2':
                DMXR2[int(pars[6:])] = pf.__dict__[pars]


    DMXList = [DMX[x] for x in DMX]
    return DMX, DMXErr, DMXR1, DMXR2



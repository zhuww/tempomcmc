'''A module that handles TOA files. 
'''
#import numpy, os
#from numpy import mean, std, var, array, double, angle
from math import *
#from fileio import *
#from SaveLoadable import MetaSaveLoader
#from decimal import *
#from numpy import float64 as __Decimal
import matplotlib as mpl 
mpl.rcParams['font.family'] = 'serif'
from decimal import Decimal as __Decimal
from decimal import getcontext
getcontext().prec = 25
secperday = 86400
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

from scipy.special import fdtr
def Ftest(chi2_1, dof_1, chi2_2, dof_2):
    """
    Ftest(chi2_1, dof_1, chi2_2, dof_2):
        Compute an F-test to see if a model with extra parameters is
        significant compared to a simpler model.  The input values are the
        (non-reduced) chi^2 values and the numbers of DOF for '1' the
        original model and '2' for the new model (with more fit params).
        The probability is computed exactly like Sherpa's F-test routine
        (in Ciao) and is also described in the Wikipedia article on the
        F-test:  http://en.wikipedia.org/wiki/F-test
        The returned value is the probability that the improvement in
        chi2 is due to chance (i.e. a low probability means that the
        new fit is quantitatively better, while a value near 1 means
        that the new model should likely be rejected).
    """
    delta_chi2 = chi2_1 - chi2_2
    delta_dof = dof_1 - dof_2
    new_redchi2 = chi2_2 / dof_2
    F = (delta_chi2 / delta_dof) / new_redchi2
    return 1.0 - fdtr(delta_dof, dof_2, F)

Observatory_list = {
        '@':   'Bary Center',
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
            self._Observatory = self.line[79].strip()[0]
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
            self._Observatory = self.line[0].strip()[0]
            self.Observatory = Observatory_list[str(self._Observatory)]
            self.info.append(self.line[2:15].strip())
            self.frequency = Decimal(self.line[15:24])
            self.TOA = Decimal(self.line[24:44])
            self.TOAsigma = Decimal(self.line[44:53])
            self.info.append(self.line[53:68].strip())
            self.DMcorr = Decimal(self.line[68:77])
        elif not line[0] ==' ' and not line[1] == ' ' and line[14] == '.' and len(line.rstrip(' ')) <= 60:
            self.format = 'ITOA'
            self.info = []
            self.line = line.ljust(59, ' ')
            self._Observatory = self.line[57:59].strip()[0]
            self.Observatory = Observatory_list[str(self._Observatory)]
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
            self._Observatory = line[4].strip()
            if self._Observatory == 'ao':
                self._Observatory = '3'
            elif self._Observatory == 'gbt':
                self._Observatory = '1'
            if self._Observatory in Observatory_list:
                self.Observatory = Observatory_list[str(self._Observatory)]
            else:
                self.Observatory = str(self._Observatory)
            #self.info.extend(line[5:])
            tags = ' '.join(line[5:])
            #self.flags = {'i':flags['i']}
            import re
            p = re.compile('-(?P<par>\w+)\s+(?P<value>\S+)', re.VERBOSE)
            for m in p.finditer(tags):
                par = str(m.group('par'))
                try:
                    value = int(m.group('value'))
                except:
                    value = str(m.group('value'))
                self.flags[par] = value
            self.info.append(tags)
            self.DMcorr = Decimal(0)


        #if self.flags.has_key('EQUAD'):
            #self.TOAsigma = sqrt(self.TOAsigma**2 + self.flags['EQUAD']**2)

    def __repr__(self):
        return '%s format TOA %s(%s us) at %s MHZ' % (self.format, self.TOA, self.TOAsigma, self.frequency)
    def __str__(self):
        TOAsigma = Decimal(self.TOAsigma)
        if self.flags.has_key('EMAX'): 
            if float(TOAsigma) > float(self.flags['EMAX']):
                return ''
        if self.flags.has_key('EMIN'): 
            if TOAsigma < Decimal(self.flags['EMIN']):
                TOAsigma = Decimal(self.flags['EMIN'])
        if self.flags.has_key('EQUAD'):
            TOAsigma = sqrt(TOAsigma**2 + self.flags['EQUAD']**2)
        if self.flags.has_key('EFAC'):
            TOAsigma = Decimal(TOAsigma) * self.flags['EFAC']
        if self.flags.has_key('to'):
            ActualTOA = self.TOA + Decimal(self.flags['to'])/secperday
        else:
            ActualTOA = self.TOA



        if self.format == 'Princeton':
            #return '%s %s%8.2f %s%9.3f%s%9.5f' % (self._Observatory, self.info[0][:13].ljust(13, ' '), self.frequency, (str(ActualTOA).ljust(20, ' '))[:20], TOAsigma, self.info[1][:15].rjust(15,' '), self.DMcorr)
            return '%s %s%8.2f %s%9.3f%s%9.5f' % (self._Observatory, self.info[0][:13].ljust(13, ' '), self.frequency, str(ActualTOA.quantize(Decimal(0.00000000000001))).rjust(20, ' '), TOAsigma, self.info[1][:15].rjust(15,' '), self.DMcorr)
        if self.format == 'Parkes':
            return ' %s%8.2f %s%8.5f%8.3f%s%s' % (self.info[0].ljust(24, ' '), self.frequency,str(ActualTOA .quantize(Decimal(0.00000000000001))).rjust(21, ' '), self.DMcorr, TOAsigma, self.info[1][:8].rjust(8, ' '), self._Observatory)
        if self.format == 'ITOA':
            try:
                info = ' '.join(self.info).strip().ljust(9, ' ')
            except TypeError:
                info = ' '.join(self.info).strip()[:9]

            if info[0] == ' ' or info[1] == ' ': info = 'IT'+info[2:]
            return '%s%s%6.3f%11.3f%10.5f  %s' % (info,str(ActualTOA .quantize(Decimal(0.0000000000001))).rjust(19, ' '), TOAsigma,self.frequency, self.DMcorr, self._Observatory)

        if self.format == 'Tempo2':
            kwpars = ''
            for key in [k for k in sorted(self.flags.keys(), reverse=True) if not k == 'EQUAD' and not k == 'JUMPflag' and not k =='EMAX' and not k =='EFAC' and not k =='EMIN']:
                kwpars += ' -%s %s ' % (key, self.flags[key])
            return '%s %s %s %s %s %s' % (self.info[0], self.frequency, str(ActualTOA .quantize(Decimal(0.00000000000001))).rjust(20, ' '), str(TOAsigma.quantize(Decimal(0.001))), self._Observatory, kwpars) 

    def tempo2fmt(self):
        TOAsigma = Decimal(self.TOAsigma)
        if self.flags.has_key('EMAX'): 
            if float(TOAsigma) > float(self.flags['EMAX']):
                return ''
        if self.flags.has_key('EMIN'): 
            if self.TOAsigma < Decimal(self.flags['EMIN']):
                TOAsigma = Decimal(self.flags['EMIN'])
        if self.flags.has_key('EQUAD'):
            TOAsigma = sqrt(TOAsigma**2 + self.flags['EQUAD']**2)
        if self.flags.has_key('EFAC'):
            TOAsigma = Decimal(TOAsigma) * self.flags['EFAC']
        if self.flags.has_key('to'):
            timeoffset = Decimal(('%.11f' % float(self.flags['to'])))/secperday
            ActualTOA = self.TOA + timeoffset
            #print '%s' % (timeoffset.quantize(Decimal(0.000000000000001)))
            #print '%s' % (self.TOA.quantize(Decimal(0.00000000000001)))
        else:
            ActualTOA = self.TOA
        if self._Observatory == '1':
            self._Observatory = 'gbt'
        elif self._Observatory == '3':
            self._Observatory = 'ao'
        file = self.file.split('/')[-1].replace(' ','_').replace('-', '_')
        file = file.ljust(10)[-15:]
        #print ActualTOA
        #fmtstr = '%s %s %s %s %s' % (file, self.frequency, str(ActualTOA.quantize(Decimal(0.00000000000001))).ljust(21, ' '), TOAsigma, self._Observatory)
        fmtstr = '%s %s %s %s %s' % (file, self.frequency, str(ActualTOA).ljust(21, ' ')[:21], TOAsigma, self._Observatory)
        kwpars = ''
        for key in [k for k in sorted(self.flags.keys(), reverse=True) if not k == 'EQUAD' and not k == 'JUMPflag' and not k =='EMAX' and not k == 'EFAC' and not k == 'EMIN' and not k == 'to']:
            kwpars += ' -%s %s ' % (key, self.flags[key])
        result = fmtstr+kwpars
        if result[14] == '.':
            result =  result[:14] + '_' + result[15:]
        return result


class TOAcommand(object):
    def __init__(self, cmd, *args):
        self.cmd = cmd
        self.args = args
    def __str__(self):
        return '%s %s' % (self.cmd, ' '.join(self.args))
    def __repr__(self):
        return '%s %s' % (self.cmd, ' '.join(self.args))

class TOAfile(object):
    """A class for read/operate TOA files. """ 
    def __init__(self, file, kws='', F0=None): 
        if kws == '':
            kws = {'JUMPflag':False} #use to be in __init()__; but I forgto why I did that.
        self.toafile = file
        self.list = []
        self.cmdlist = []
        self.PHASEJUMPS = {}
        self.EQUADvalues = {}
        self.start = Decimal('1000000')
        self.end = Decimal('0')
        self.F0 = F0
        with open(file,'rt') as f:
            text = f.read()
            lines = text.split('\n')
            toagroup = []
            SKIPflag = False
            PHASEJUMPFlag = False
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
                    self.cmdlist.extend(tf.cmdlist)
                elif len(l) < 25:
                    if len(toagroup) > 0:self.list.append(toagroup)
                    command = TOAcommand(*l.split())
                    self.cmdlist.append(command)
                    toagroup = []
                    if command.cmd == 'INFO':
                        INFO = command.args[0]
                        kws.update({'i': INFO})
                        if PHASEJUMPFlag:
                            self.PHASEJUMPS[INFO] = PHASE
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
                        if kws.has_key('padd'):
                            NewPhase = kws['padd'] + PHASE
                            if NewPhase == 0.:
                                del kws['padd']
                                PHASEJUMPFlag = False
                            else:
                                kws['padd'] = NewPhase
                        else:
                            kws['padd'] = PHASE
                            PHASEJUMPFlag = True
                    elif command.cmd == 'EQUAD':
                        EQUAD = Decimal(command.args[0])
                        kws.update({'EQUAD':EQUAD})
                        if not kws.has_key('f'):
                            key = 1
                            self.EQUADvalues[str(key)]=EQUAD
                        else:
                            key = int(kws['f']) + 1
                            self.EQUADvalues[str(key)]=EQUAD
                        kws.update({'f':str(key)})
                    elif command.cmd == 'EMAX':
                        self.EMAX = Decimal(command.args[0])
                        kws['EMAX'] = self.EMAX
                        self.list.append(command)#ignore
                    elif command.cmd == 'EMIN':
                        self.EMIN = Decimal(command.args[0])
                        kws['EMIN'] = self.EMIN
                        self.list.append(command)#ignore
                    elif command.cmd == 'EFAC':
                        self.EFAC = Decimal(command.args[0])
                        self.list.append(command)#ignore
                        kws['EFAC'] = self.EFAC
                    elif command.cmd == 'SKIP':
                        SKIPflag = True
                    elif command.cmd == 'NOSKIP':
                        SKIPflag = False
                    else:
                        self.list.append(command)#ignore
                else:
                    if len(self.cmdlist) == 0 or not self.cmdlist[-1] == 'TOAlist':
                        self.cmdlist.append('TOAlist')
                    if not SKIPflag:
                        if not file.find(' ') == -1:
                            file = ''
                        try:
                            t= TOA(l, file = file, **kws)
                        except TypeError:
                            print 'Probamatic line: ', l
                            print 'from file %s' % file
                        if t.flags.has_key('EMAX'): 
                            if float(t.TOAsigma) > float(t.flags['EMAX']):
                                continue
                        #if t.flags.has_key('EFAC') and t.flags.has_key('EMAX'):
                            #if t.TOAsigma * t.flags['EFAC'] > t.flags['EMAX']:
                                #continue 

                        toagroup.append(t)
                        if t.TOA < self.start:self.start = t.TOA
                        if t.TOA > self.end: self.end = t.TOA
                    else:
                        pass # skipped
            if len(toagroup) > 0:self.list.append(toagroup)
        toalist = []
        for toa in self.list:
            if not type(toa) == TOAcommand:
                toalist.extend(toa)
        self.toalist = toalist
        self.issubgroup = False
        self.inspect()

    def inspect(self):
        commandgrps = {}
        subcmdgrp = []
        info = 'Untaged'
        for c in self.cmdlist:
            if type(c) == TOAcommand:
                if not c.cmd.find('SKIP') == -1:continue
                if c.cmd in ['EQUAD', 'EMAX', 'EFAC', 'EMIN']:
                #if c.cmd in ['EQUAD', 'EMAX']:
                    #subcmdgrp.append(c)
                    pass
                if c.cmd == 'INFO':info = c.args[0]
                if c.cmd == 'MODE':self.MODE = c.args[0]
            else:
                #if not len(subcmdgrp) == 0:
                if not commandgrps.has_key(info):
                    commandgrps[info] = subcmdgrp
                else:
                    commandgrps[info].extend(subcmdgrp)
                subcmdgrp = []
        self.commandgrps = commandgrps
        toagrp = {}
        jumpgrp = {}
        phasegroups = {}
        othergroups = {}
        toalist = self.toalist
        for i in range(len(toalist)):
            toa = toalist[i]
            #if not toa.flags.has_key('i'):
                #info = 'Untaged'
            #else:
                #info = toa.flags['i']
            #if not toa.flags.has_key('jump'):
                #jump = 0
            #else:
                #jump = toa.flags['jump']
            info = 'Untaged'
            jump = 0
            for key in toa.flags:
                if key == 'i':
                    info = toa.flags['i']
                elif key == 'jump':
                    jump = toa.flags['jump']
                elif key == 'padd':
                    phase = float(toa.flags['padd'])

                elif toa.flags.has_key('padd'):
                    phase = float(toa.flags['padd'])
                    if phasegroups.has_key(info):
                        if not phasegroups[info] == phase:raise "different phase jump for the same info group %s" % info
                    else:
                        phasegroups[info] = phase
                else:
                    keyword = '-'+key+' '+str(toa.flags[key])
                    if othergroups.has_key(keyword):
                        othergroups[keyword].append(i)
                    else:
                        othergroups[keyword] = [i]
            if not toagrp.has_key(info):
                toagrp[info] = [i]
            else:
                toagrp[info].append(i)
            if not jumpgrp.has_key(jump):
                jumpgrp[jump] = [i]
            else:
                jumpgrp[jump].append(i)
        self.groups = toagrp
        self.jumpgroups = jumpgrp
        self.phasegroups = phasegroups
        self.othergroups = othergroups
        #print jumpgrp

        self.matchdict = {}
        self.matchnotag = {}
        self.firstgrp = None
        self.grouporder = []
        for grp in self.groups.keys():#now go through all groups and see if they contain the jump gropups
            self.matchdict[grp] = []
            if self.jumpgroups.has_key(0):
                if set(self.jumpgroups[0]) <= set(self.groups[grp]):
                    self.matchdict[grp].append(0)
                    self.firstgrp = grp
                    self.grouporder.append(grp)
        if self.firstgrp == None:
            self.matchdict['notag'] = [0]
            self.firstgrp = 'notag'
            self.grouporder.append('notag')

        for jpgrp in [k for k in self.jumpgroups.keys() if not k == 0]:
            tagged = False
            for grp in self.groups.keys():
                if set(self.jumpgroups[jpgrp]) <= set(self.groups[grp]) and tagged == False:
                    self.matchdict[grp].append(jpgrp)
                    tagged = True
                    if not self.grouporder[-1] == grp and not grp in self.grouporder:
                        self.grouporder.append(grp)
            if tagged == False:
                if not self.matchdict.has_key('notag'):
                    self.matchdict['notag'] = [jpgrp]
                else:
                    self.matchdict['notag'].append(jpgrp)
                if not self.grouporder[-1] == 'notag' and not 'notag' in self.grouporder:
                    self.grouporder.append('notag')
        if self.matchdict.has_key('notag') and any([len(self.matchdict[k])==0 for k in self.matchdict.keys()]):
            self.matchnotag = {} #See if jump group contains multiple info tags
            for jpgrp in self.matchdict['notag']:
                for grp in [k for k in self.matchdict.keys() if len(self.matchdict[k]) == 0 ]:
                    #print grp, self.matchdict[grp],len(self.matchdict[grp])
                    if self.jumpgroups.has_key(0) and set(self.groups[grp]) < set(self.jumpgroups[jpgrp]):
                        if not self.matchnotag.has_key(jpgrp):
                            self.matchnotag[jpgrp] = [grp]
                        else:
                            self.matchnotag[jpgrp].append(grp)
            alluntagged =set()
            #print self.toafile
            #print self.matchnotag
            #print self.matchdict
            #print self.jumpgroups.keys()
            for jpgrp in self.matchdict['notag']:
                #print jpgrp, self.matchnotag[jpgrp]
                for s in [set(self.groups[grp]) for grp in self.matchnotag[jpgrp]]:
                    alluntagged |= s
            alluntaggedjp =set()
            for s in [set(self.jumpgroups[jpgrp]) for jpgrp in self.matchdict['notag']]:
                alluntaggedjp |= s
            if not  alluntagged == alluntaggedjp:
                print 'self.matchnotag: ', self.matchnotag
                print 'self.matchdic: ', self.matchdict
                #print [grp for grp in self.matchnotag[jpgrp] for jpgrp in [0,1]]
                #print [grp for grp in self.matchnotag[jpgrp] for jpgrp in self.matchdict['notag']]
                #print alluntagged #, alluntaggedjp
                #print alluntaggedjp
                print 'Warning: some toas does not belong to a jump group:' ,
                self.groups['UntagedJump'] = list(alluntaggedjp - alluntagged)
                print self.groups['UntagedJump']

            self.jumpdict = self.matchdict.copy()
            for jumpgrp in self.matchnotag.keys():
                for grp in self.matchnotag[jumpgrp]:
                    self.jumpdict[grp] = [jumpgrp]
            self.jumpdict.pop('notag')
        if not self.toalist[0].format == 'Tempo2':
            self.format = 'Tempo1'
        else:
            self.format = 'Tempo2'


    def subgroup(self, groups=None):
        from copy import deepcopy
        new = deepcopy(self)
        if groups == None:
            return new
        elif not type(groups) == list:
            raise "Expecting groups to be list by get %s" % type(groups)
        new.start = Decimal('1000000')
        new.end = Decimal('0')
        new.list = []
        nextjump = 0
        newjumpgrp = {}
        subcmdgrp = []
        subtoalist = []
        EQUADgroups = []
        for c in self.list:
            if type(c) == TOAcommand:
                if not c.cmd.find('SKIP') == -1:continue
                if c.cmd in ['EQUAD', 'EMAX', 'EFAC', 'EMIN']:
                    subcmdgrp.append(c)
                if c.cmd == 'INFO':
                    info = c.args[0]
                    subcmdgrp.append(c)
                if not len(subtoalist) == 0:
                    new.list.append(subtoalist)
                    subtoalist = []
            else:
                if not len(subcmdgrp) == 0:
                    if info in groups:
                        new.list.extend(subcmdgrp)
                    subcmdgrp = []
                for t in c:
                    if t.flags.has_key('i') and t.flags['i'] in groups:
                        if t.flags.has_key('jump'):
                            i = int(t.flags['jump'])
                        else:
                            i = 0
                        if newjumpgrp.has_key(i):
                            if newjumpgrp[i] == 0:
                                try:
                                    t.flags.pop('jump')
                                except:pass
                            else:
                                t.flags['jump'] = newjumpgrp[i]
                        else:
                            newjumpgrp[i] = nextjump
                            nextjump += 1
                            t.flags['jump'] = newjumpgrp[i]
                        subtoalist.append(t)
                        if t.TOA < new.start:new.start = t.TOA
                        if t.TOA > new.end: new.end = t.TOA
                        if t.flags.has_key('f') and not t.flags['f'] in EQUADgroups:
                            EQUADgroups.append(t.flags['f'])
                    else:
                        #print t.flags
                        pass
        if not len(subtoalist) == 0:
            new.list.append(subtoalist)
            subtoalist = []
        toalist = []
        for toa in new.list:
            if not type(toa) == TOAcommand:
                toalist.extend(toa)
        new.toalist = toalist
        for key in new.EQUADvalues.keys():
            if not key in EQUADgroups:
                new.EQUADvalues.pop(key)
        new.inspect()
        new.newjumpgrp = newjumpgrp
        new.newjumpgroup = {}
        for key in newjumpgrp:
            new.newjumpgroup[newjumpgrp[key]] = key
        if new.toalist[0].__dict__.has_key('npulse'):
            new.npulse = []
            for i in range(len(new.toalist)):
                t = new.toalist[i]
                new.npulse.append(t.npulse)
        new.issubgroup = True
        return new


    def __str__(self):
        #text = ''
        #for x in self.list:
            #if type(x) == TOAcommand:
                #text += str(x)+'\n'
            #else:
                #for t in x:
                    #text += str(t)+'\n'
        #return text
        return  '\n'.join([str(x) for x in self.list])


    def tempo1fmt(self, format = 'Princeton'):
        """convert TOA file to tempo1 format, use the format keyword to define which version (Princeton, Parkes or ITOA) """
        for toa in self.toalist:
            toa.format = format
        text = ''
        if self.__dict__.has_key('MODE'):
            text += 'MODE %s\n' % self.MODE
        phaseoffset = 0. 
        for grp in self.grouporder:
            if not grp == 'notag':
                if grp in self.phasegroups.keys():
                    if not phaseoffset == self.phasegroups[grp]:
                        text += 'PHASE %s\n' % (self.phasegroups[grp] - phaseoffset)
                        phaseoffset = self.phasegroups[grp]
                else:
                    if not phaseoffset == 0.:
                        text += 'PHASE %s\n' % (0. - phaseoffset)
                        phaseoffset = 0.

                text += 'INFO %s\n' % (grp)
                for cmd in self.commandgrps[grp]:
                    text += '%s\n' % cmd
                for jpgrp in self.matchdict[grp]:
                    if not jpgrp == 0:
                        text += 'JUMP\n'
                    for i in self.jumpgroups[jpgrp]:
                        toa = self.toalist[i]
                        oldfmt = toa.format
                        toa.format = format
                        text += '%s\n' % str(toa)
                        toa.format = oldfmt
                    if not jpgrp == 0:
                        text += 'JUMP\n'
            else:
                for jpgrp in self.matchdict['notag']:
                    if not jpgrp == 0:
                        text += 'JUMP\n'
                    for grp in self.matchnotag[jpgrp]:
                        if grp in self.phasegroups.keys():
                            if not phaseoffset == self.phasegroups[grp]:
                                text += 'PHASE %s\n' % (self.phasegroups[grp] - phaseoffset)
                                phaseoffset = self.phasegroups[grp]
                        else:
                            if not phaseoffset == 0.:
                                text += 'PHASE %s\n' % (0. - phaseoffset)
                                phaseoffset = 0.
                        text += 'INFO %s\n' % (grp)
                        for cmd in self.commandgrps[grp]:
                            text += '%s\n' % cmd
                        for i in self.groups[grp]:
                            toa = self.toalist[i]
                            oldfmt = toa.format
                            toa.format = format
                            text += '%s\n' % str(toa)
                            toa.format = oldfmt
                    if not jpgrp == 0:
                        text += 'JUMP\n'

        return text



    def tempo2fmt(self, tempo1use=False):
        """convert TOA file to tempo2 format"""
        fmtstr = """
"""
        fmtstr += '\n'
        hasphasejump = False
        currentphasejump = 0
        accumulatephasejump = 0
        jumpnumber = 0
        hasjump = False
        '''*** try to temporarily put PHASE/jumps in tempo2 format for tempo1 use'''
        for toas in self.list:
            if isinstance(toas, (list,tuple)):
                if tempo1use:
                    if 'padd' in toas[0].flags and not hasphasejump:
                        currentphasejump = toas[0].flags['padd']
                        accumulatephasejump += toas[0].flags['padd']
                        fmtstr += 'PHASE %s\n' % (currentphasejump)
                        hasphasejump = True
                    elif 'padd' in toas[0].flags and hasphasejump:
                        currentphasejump = toas[0].flags['padd'] - accumulatephasejump
                        if not currentphasejump == 0.:
                            fmtstr += 'PHASE %s\n' % (currentphasejump)
                    elif not 'padd' in toas[0].flags and hasphasejump:
                        currentphasejump = 0 - accumulatephasejump
                        fmtstr += 'PHASE %s\n' % (currentphasejump)
                        hasphasejump = False

                    if 'jump' in toas[0].flags and jumpnumber == 0:
                        jumpnumber = toas[0].flags['jump']
                        fmtstr += 'JUMP\n'
                    elif 'jump' in toas[0].flags and not toas[0].flags['jump'] == jumpnumber:
                        jumpnumber = toas[0].flags['jump']
                        fmtstr += 'JUMP\n'
                        fmtstr += 'JUMP\n'

                '''*** try to temporarily put PHASE/jumps in tempo2 format for tempo1 use'''

                fmtstr += '\n'.join([t.tempo2fmt() for t in toas])
                fmtstr += '\n'
            elif isinstance(toas, TOAcommand) and not toas.cmd in ['EMAX', 'EFAC', 'EQUAD', 'MODE', 'EMIN']:
                fmtstr += '%s\n' % (toas)


        Extrastr = 'FORMAT 1\nMODE 1\n'
        #for key in self.EQUADvalues.keys():
            #Extrastr += 'T2EQUAD -f %s %s\n' % (key, self.EQUADvalues[key])
        return Extrastr + fmtstr


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
    """
    tempofit(parfile, toafile=None, pulsefile=None):
    use tempo to fit the parfile and toafile, specify pulsefile to use pulse number.
    """
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
        #print parfile, toafile
        pass
    try:
        dof = int(line.split('/')[2].strip(' ').split()[0])
    except:
        #print line
        dof = 13722
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
        'FD1  ':'FD1',
        'FD2  ':'FD2',
        'FD3  ':'FD3',
        'FD4  ':'FD4',
        'FD*  ':'FD5',
        'DM001':'DM001',
        'DM002':'DM002',
        'eps1 ':'EPS1',
        'eps2 ':'EPS2',
        }
for i in range(500):
    key = 'DX' + str(i).rjust(3,'0')
    value = 'DMX_' + str(i).rjust(4,'0')
    paramap[key] = value
for i in range(500):
    key = 'D1' + str(i).rjust(3,'0')
    value = 'DMX1_' + str(i).rjust(4,'0')
    paramap[key] = value
for i in range(10,1000):
    key = ' O' + str(i).rjust(3,'0')
    value = 'JUMP_' + str(i)
    paramap[key] = value

class PARfile(object):
    #__metaclass__ = MetaSaveLoader
    def __init__(self, file):
        """A class object to contain the APIs for extracting ephemeris information from a .par file."""
        if os.access(file, os.R_OK):
            self.parfile = file
        else: 
            print os.getcwd() + '/' + file
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
        self.jumps = {0:Decimal(0)}
        self.TempoVersion = 1
        for lines in self.file.readlines():
            items = lines.split()
            if len(items) == 0:pass
            elif items[0] == '#':pass
            elif items[0][0] == '#':pass
            elif items[0]== 'C':pass
            elif items[0].upper() == 'JUMP' and items[1].startswith('-'):
                jumptag = '%s %s %s' % tuple(items[:3])
                self.manifest.append(jumptag)
                if len(items) > 3:
                    self.__dict__[jumptag] = [Decimal(items[3]), Decimal(0)]
                    if len(items) > 4:
                        if items[4] == '0':
                            self.parameters[jumptag] = '0'
                        else:
                            self.parameters[jumptag] = '1'
                    if len(items) > 5:
                        self.__dict__[jumptag][1] = Decimal(items[5])
                else:
                    self.__dict__[jumptag] = (Decimal(0), Decimal(0))
                    self.parameters[jumptag] = '1'
            elif len(items) == 2 and not items[0] in ['SINI', 'M2', 'XDOT']: 
                self.__dict__[items[0]] = floatify(items[1])
                self.manifest.append(items[0])
            elif items[0] in ['START', 'FINISH']:
                if len(items)>2 and items[2] == '1':
                    self.__dict__[items[0]] = items[1] + ' ' + '1'
                else:
                    self.__dict__[items[0]] = value
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
                if not items[0].find('JUMP') == -1:
                    if not items[0].find('_') == -1:
                        try:
                            self.jumps[int(items[0].split('_')[-1])] = value
                        except:
                            self.jumps[items[0].split('_')[-1]] = value
                    else:
                        try:
                            self.jumps[int(items[2])] = floatify(items[3])
                        except:
                            self.jumps[str(items[2])] = floatify(items[3])
            elif len(items) >= 5 and items[0] == 'JUMP':
                name = ' '.join(items[:3])
                self.manifest.append(name)
                value = floatify(items[3])
                fitflag = items[4]
                error = ''
                self.__dict__[name] = [value, error]
                self.parameters.update({name:fitflag})
                try:
                    if items[1] == '-jump':
                        self.jumps[int(items[2])] = value
                    else:
                        self.jumps[' '.join(items[1:3])] = value
                except:
                    self.jumps[' '.join(items[1:3])] = value
                self.TempoVersion = 2
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
            self.psrname = self.PSRJ
        elif self.__dict__.has_key('PSRB'): 
            self.psrname = self.PSRB
        else:
            raise 'Cant Find PSR name'
        self.file.close()
        if self.__dict__.has_key('BINARY') and self.BINARY == 'T2':self.TempoVersion = 2
    def write(self, fn=None):
        if fn == None:
            fn = self.parfile
        with open(fn, 'w') as f:
            text = ''
            for item in self.manifest:
                if item in self.parameters.keys() :
                    if item.startswith('DMX'):
                        text += '%s\t%s\t%s\t%s\n' % (item, self.__dict__[item][0], self.parameters[item] ,self.__dict__[item][1])#).replace('E+', 'D+').replace('E-','D-')
                    else:
                        text += '%s\t%s\t%s\t%s\n' % (item, self.__dict__[item][0], self.parameters[item] ,self.__dict__[item][1])#).replace('E+', 'D+').replace('E-','D-')
                else:
                    text += ('%s\t%s\n' % (item, self.__dict__[item])).replace('D+00','')
            f.write(text)
    def freezeall(self, key=None):
        import re
        if key == None:
            for p in self.parameters:
                self.parameters[p] = '0'
        else:
            for p in self.parameters:
                if not re.match(key, p) == None:
                    #print 'freezed parameter %s' % p
                    self.parameters[p] = '0'
    def thawall(self, key = None):
        import re
        for p in self.parameters:
            if key == None:
                self.parameters[p] = '1'
            else:
                if not re.match(key, p) == None:
                    self.parameters[p] = '1'

    def matrix(self, timfile, TempoVersion=1):
        if self.BINARY == 'DDS':
            paramap['    s'] = 'SHAPMAX'
        if TempoVersion == 1:
            from numpy.core.records import fromfile
            os.system('tempo -j -f %s %s > tmplog' % (self.parfile, timfile))
            bpf = PARfile(self.psrname + '.par')
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
                if not paramj == 'FD*  ':
                    try:
                        param = paramap[paramj]
                    except KeyError:
                        param = paramj
                    self.parlist.append(param)
                #try:
                if param == 'T0':
                    if 'BINARY' in bpf.__dict__ and str(bpf.BINARY) == 'ELL1':
                        param = 'TASC'
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
                self.Nfac = int(sqrt(m))+1
                self.par = [0]*m
                #print self.covariance[:2,:2]
        else:
            print 'TempoVersion: ', TempoVersion
            raise TempoVersion
    def randomnew(self, stepsize=1.):
        new = deepcopy(self)
        err = mvn(self.par, self.covariance)
        for i in range(len(self.parlist)):
            p = self.parlist[i]
            if p == 'RAJ' or p == 'DECJ':
                dd,mm,ss = new.__dict__[p][0].split(':')
                #ss = str(Decimal(ss) + Decimal(str(err[i]*float(str(new.__dict__[p][1]))/self.Nfac)))
                #ss = str(Decimal(ss) + Decimal(repr(err[i]*self.err[i]/self.Nfac)))
                ss = str(Decimal(ss) + Decimal(repr(err[i]*self.err[i]/self.Nfac*stepsize)))
                new.__dict__[p][0] = ':'.join([dd,mm,ss])
            elif p == 'SINI':
                pass
            #elif p == 'PAASCNODE':
                #new.__dict__[p][0] = new.__dict__[p][0] + Decimal(np.random.rand()*10/self.Nfac*stepsize)
            else:
                #new.__dict__[p][0] = new.__dict__[p][0] + Decimal(str(err[i]*float(str(new.__dict__[p][1]))/self.Nfac))
                new.__dict__[p][0] = new.__dict__[p][0] + Decimal(repr(err[i]*self.err[i]/self.Nfac*stepsize))

        return new
    def convert(self, TempoVersion):
        """Convert tempo1 parfile to tempo2 version, or vise versa """
        if self.TempoVersion == TempoVersion:pass
        else:
            if TempoVersion == 1:
                self.TempoVersion = 1
                for i in range(len(self.manifest)):
                    par = self.manifest[i]
                    if not par.find('JUMP') == -1:
                        j = par.split()[-1]
                        newjumppar = 'JUMP_' + j
                        self.manifest[i] = newjumppar
                        self.__dict__[newjumppar] = self.__dict__[par]
                        del self.__dict__[par]
                        self.parameters[newjumppar] = self.parameters[par]
                        del self.parameters[par]
                if not self.CLK.find('UTC') == -1:
                    self.CLK = 'TT(BIPM)'
                if self.BINARY == 'T2':
                    self.BINARY = 'DD'
                    self.E = self.ECC
                    del self.__dict__['ECC']
                    self.parameters['E'] = self.parameters['ECC']
                    del self.parameters['ECC']
                    i = self.manifest.index('ECC')
                    self.manifest[i] = 'E'
                    self.PAASCNODE = self.KOM[0] 
                    del self.__dict__['KOM']
                    i = self.manifest.index('KOM')
                    self.manifest[i] = 'PAASCNODE'
                    del self.parameters['KOM']
                    self.SINI = [sin(float(self.KIN[0])/180.*pi), abs(cos(float(self.KIN[0])/180.*pi)*float(self.KIN[1])/180.*pi)]
                    i = self.manifest.index('KIN')
                    #self.manifest[i] = 'SINI'
                    self.manifest.pop(i)
                    inc = float(self.KIN[0])/180.*pi
                    self.parameters['SINI'] = 1
                    del self.parameters['KIN']
                    del self.__dict__['KIN']
                    fac = 0.0001536
                    mu = sqrt(float(self.PMRA[0])**2 + float(self.PMDEC[0])**2)
                    thetamu = 90. - atan(float(self.PMDEC[0])/float(self.PMRA[0]))*180./pi
                    Omega = float(self.PAASCNODE)
                    x = float(self.A1[0])
                    #print 'mu, thetamu, x', mu, thetamu, x
                    xdot = fac * x * mu * (cos(inc)/sin(inc)) * sin((thetamu-Omega)*pi/180.)
                    self.__dict__['XDOT'] = (xdot, 0.1*xdot)
                    self.manifest.append('XDOT')
                    self.parameters['XDOT'] = '1'
                    for par in ['EPHVER', 'MODE', 'UNITS', 'TIMEEPH', 'DILATEFREQ', 'PLANET_SHAPIRO', 'T2CMETHOD', 'CORRECT_TROPOSPHERE', 'CHI2R', 'TRES', 'NTOA', 'POSEPOCH']:
                        try:
                            self.manifest.pop(self.manifest.index(par))
                        except:pass
                    if 'DMX_0001' in self.manifest:
                        for par in self.manifest:
                            if not par.find('DMX_') == -1:
                                par1 = 'DMX1_' + par.split('_')[-1]
                                self.__dict__[par1] = (Decimal(0), Decimal(0))
                                self.parameters[par1] = '0'
                                self.manifest.append(par1)
                    DMXlist = ['DMX']
                    self.__dict__['DMX'] = '0.10000000D+02'
                    JUMPlist = []
                    manifest = self.manifest[:]
                    for par in self.manifest:
                        if not par.find('DMX_') == -1:
                            i = par.split('_')[-1]
                            DMXlist.append('DMX_'+i)
                            DMXlist.append('DMX1_'+i)
                            DMXlist.append('DMXR1_'+i)
                            DMXlist.append('DMXR2_'+i)
                            #DMXlist.append('DMXF1_'+i)
                            #DMXlist.append('DMXF2_'+i)
                            #self.__dict__['DMXF1_'+i] = Decimal('700.')
                            #self.__dict__['DMXF2_'+i] = Decimal('3100.')
                            manifest.pop(manifest.index(par))
                        elif not par.find('DMX') == -1:
                            manifest.pop(manifest.index(par))
                        elif not par.find('JUMP') == -1:
                            manifest.pop(manifest.index(par))
                            JUMPlist.append(par)
                        else:
                            pass
                    self.manifest = manifest
                    JUMPlist.sort(key = lambda x:int(x.split('_')[-1]))
                    self.manifest.extend(DMXlist)
                    self.manifest.extend(JUMPlist)
                    #self.manifest = list(set(self.manifest))
                    
            else:
                self.TempoVersion = 2
                for i in range(len(self.manifest)):
                    par = self.manifest[i]
                    if not par.find('JUMP') == -1:
                        if not par.find(' -') == -1:
                            newjumppar = par
                        else:
                            j = par.split('_')[-1]
                            newjumppar = 'JUMP -jump ' + j
                        self.manifest[i] = newjumppar
                        if not newjumppar == par:
                            self.__dict__[newjumppar] = self.__dict__[par]
                            del self.__dict__[par]
                            self.parameters[newjumppar] = self.parameters[par]
                            del self.parameters[par]
                if self.BINARY == 'DD' and self.__dict__.has_key('PAASCNODE'):
                    self.BINARY = 'T2'
                    self.ECC = self.E
                    del self.__dict__['E']
                    self.parameters['ECC'] = self.parameters['E']
                    del self.parameters['E']
                    i = self.manifest.index('E')
                    self.manifest[i] = 'ECC'
                    self.KOM = (self.PAASCNODE, self.PAASCNODE/10)
                    del self.__dict__['PAASCNODE']
                    i = self.manifest.index('PAASCNODE')
                    self.manifest[i] = 'KOM'
                    #del self.parameters['PAASCNODE']
                    self.parameters['KOM'] = '1'
                    #self.SINI = [sin(float(self.KIN[0])/180.*pi), abs(cos(float(self.KIN[0])/180.*pi)*float(self.KIN[1])/180.*pi)]
                    i = self.manifest.index('SINI')
                    self.manifest[i] = 'KIN'
                    inc = asin(float(self.SINI[0]))*180./pi
                    self.parameters['KIN'] = 1
                    self.__dict__['KIN'] = [inc, float(self.__dict__['SINI'][1])/cos(inc/180.*pi)*180./pi]
                    del self.parameters['SINI']
                    self.SINI = 'KIN'
                    del self.__dict__['XDOT']
                    self.manifest[self.manifest.index('XDOT')] = 'SINI'
                    #fac = 0.0001536
                    #mu = sqrt(float(self.PMRA[0])**2 + float(self.PMDEC[0])**2)
                    #thetamu = 90. - atan(float(self.PMDEC[0])/float(self.PMRA[0]))*180./pi
                    #Omega = float(self.PAASCNODE)
                    #x = float(self.A1[0])
                    #print 'mu, thetamu, x', mu, thetamu, x
                    #xdot = fac * x * mu * (cos(inc)/sin(inc)) * sin((thetamu-Omega)*pi/180.)
                    #self.__dict__['XDOT'] = (xdot, 0.1*xdot)
                    #self.manifest.append('XDOT')
                    #self.parameters['XDOT'] = '1'







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
    if type(pf.__dict__['DM']) in (list, tuple):
        DM = pf.__dict__['DM'][0]
    else:
        DM = pf.__dict__['DM']
    for pars in DMXpars:
        if pars.find('DMXR') == -1 and len(pars)>=4:
            if not pars[3] in 'EF':
                if not pars[3] == '1':
                    try:
                        DMX[int(pars[4:])] = pf.__dict__[pars][0] + DM
                        DMXErr[int(pars[4:])] = pf.__dict__[pars][1] 
                    except:
                        #print pars , pf.__dict__[pars]
                        #raise stop
                        DMX[int(pars[4:])] = pf.__dict__[pars] + DM
                        DMXErr[int(pars[4:])] = 0
                else:
                    DMX1[int(pars[5:])] = pf.__dict__[pars][0]
        elif not pars.find('DMXR') == -1 and len(pars)>=4:
            if pars[4] == '1':
                DMXR1[int(pars[6:])] = pf.__dict__[pars]
            elif pars[4] == '2':
                DMXR2[int(pars[6:])] = pf.__dict__[pars]


    #DMXList = [DMX[x] for x in DMX]
    return DMX, DMXErr, DMXR1, DMXR2


import tempfile
class model(PARfile):
    """ An extension of the PARfile object, which stores the methods and results of a tempo 1 or 2 fitting. """

    def __init__(self, parfile):
        super(self.__class__, self).__init__(parfile)
        self.NITS = '1'
        tmpparfile = tempfile.NamedTemporaryFile(mode='w', suffix=".par", delete=False).name
        self.write(tmpparfile)
        self.parfile = tmpparfile
        if self.__dict__.has_key('DMX') or self.__dict__.has_key('DMX_0001'):
            self.dmxlist = getDMX(self)

    def tempofit(self, toafile, pulsefile=None):
        """
        model.tempofit(toafile, pulsefile=None):
        use tempo to fit the parfile and toafile, specify pulsefile to use pulse number.
        """
        from numpy import mean, array, fromfile, savetxt
        self.START = toafile.start
        self.FINISH = toafile.end
        self.toafile = toafile
        if self.TempoVersion == 2:
            self.convert(1)
            self.write()
        tmppulsefile = tempfile.NamedTemporaryFile(mode='w', suffix='.pls', delete = False).name
        #if toafile.issubgroup or not toafile.format == 'Tempo1':
        if toafile.issubgroup:
            tmptimfile = tempfile.NamedTemporaryFile(mode='w', suffix=".tim", delete=False).name
            tf = open(tmptimfile, 'w')
            tf.write(toafile.tempo1fmt())
            tf.close()
            toafile.toafile = tmptimfile
        if toafile.issubgroup:
            jumpdict = {}
            firstgrp = toafile.newjumpgroup[0]
            if firstgrp == 0:
                self.JUMP_0 = (Decimal(0), Decimal(0))
            self.jumps = {}
            for newjp in toafile.newjumpgroup.keys():
                if not newjp == 0:
                    oldjp = toafile.newjumpgroup[newjp]
                    try:
                        jumpdict['JUMP_%d' % newjp] = (self.__dict__['JUMP_%d' % oldjp][0] - self.__dict__['JUMP_%d' % firstgrp][0], self.__dict__['JUMP_%d' % oldjp][1] + self.__dict__['JUMP_%d' % firstgrp][1] )
                    except:
                        print oldjp, newjp
                        print jumpdict
                        print self.__dict__['JUMP_%d' % oldjp]
                    self.jumps[newjp] = jumpdict['JUMP_%d' % newjp][0]
                    
            if firstgrp == 0:
                del self.__dict__['JUMP_0'] 
            for key in self.__dict__.keys():
                if not key.find('JUMP_') == -1:
                    del self.__dict__[key]
                    self.manifest.pop(self.manifest.index(key))
                    del self.parameters[key] #clean up old jumps
            for key in jumpdict.keys():
                self.__dict__[key] = jumpdict[key]
                self.manifest.append(key)
                self.parameters[key] = '1' #default to be fitted

        if self.__dict__.has_key('DMX') or self.__dict__.has_key('DMX_0001'):
            DMX, DMXErr, DMXR1, DMXR2 = self.dmxlist   
            for key in DMXR1:
                if DMXR2[key] < self.START or DMXR1[key] > self.FINISH:
                    self.parameters['DMX_' + str(key).rjust(4,'0')] = '0'
                else:
                    pass
                #self.parameters['DMX_' + str(key).rjust(4,'0')] = '1'

        self.groups = toafile.groups.copy()
        self.jumpgroups = toafile.jumpgroups.copy()
        self.write()
        #if not set(self.jumps.keys()) <= set(toafile.jumpgroups.keys()):
            #print "Jump group flag mismatch!"
        if pulsefile == None:
            if toafile.__dict__.has_key('npulse'):
                savetxt(tmppulsefile, toafile.npulse, fmt='%.0f')
                line = getoutput("tempo  -f %s %s -ni %s -j| grep 'Chisqr/nfree'" % (self.parfile, toafile.toafile, tmppulsefile))
                #print 'tempo  -f %s %s -ni %s -j' % (self.parfile, toafile.toafile, tmppulsefile)
            else:
                line = getoutput("tempo -f %s %s -j -no %s| grep 'Chisqr/nfree'" % (self.parfile, toafile.toafile, tmppulsefile))
                #print tmppulsefile
        else:
            line = getoutput("tempo  -f %s %s -ni %s -j| grep 'Chisqr/nfree'" % (self.parfile, toafile.toafile, pulsefile))
            tmppulsefile = pulsefile
        self.line = line
        try:
            chisq = line.split('/')[1].split(':')[1]
        except:
            pass
        try:
            dof = int(line.split('/')[2].strip(' ').split()[0])
        except:
            dof = len(toafile.list) - len([ x for x in self.parameters])
        try:
            chisq = float(chisq)
        except:
            chisq = 9999999.
        self.chisq = chisq 
        self.dof = dof
        #try:
        self.newpar = PARfile(self.psrname + '.par')
        #except IndexError:
            #print self.parfile
            #print self.toafile.toafile

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

        self.toa = data[:]['toa']
        self.freq = data[:]['rf_bary']
        self.res = data[:]['res_sec']*1.e6
        self.err = data[:]['err_us']
        self.phase = data[:]['res_phase']
        self.ophase = data[:]['ophase']
        self.prefitres= data[:]['prefit_sec']*1.e6
        self.prefitphase= self.prefitres*self.phase/self.res
        self.weight = data[:]['weight']
        #from fileio import readcol
        npulse= np.genfromtxt(tmppulsefile, dtype='int')
        #npulse= readcol(tmppulsefile, 0)
        self.toafile.npulse = npulse
        if self.chisq/self.dof < 2.: #if the fit is ok, then log the pulse number
            for i in range(len(self.toafile.toalist)):
                t = self.toafile.toalist[i]
                try:
                    t.npulse = npulse[i]
                except IndexError:
                    print i, len(npulse)
            self.toafile.pulsefile = tmppulsefile

        from datatools.MJD import MJD_to_datetime
        from datetime import timedelta
        bat = self.toa
        toalist = self.toafile.toalist
        for i in range(len(bat)):
            toa = toalist[i]
            if MJD_to_datetime(float(toa.TOA)) - MJD_to_datetime(bat[i])  > timedelta(seconds = 600):
                for m in range(3,1,-1):
                    print i-m, 'TOA', toalist[i-m].TOA, bat[i-m], (MJD_to_datetime(float(toalist[i-m].TOA)) - MJD_to_datetime(bat[i-m])).seconds, toalist[i-m].TOAsigma, self.err[i-m], toalist[i-m].flags['i'], toalist[i-m].flags['EMAX']
                print '>', i, 'TOA', toa.TOA, bat[i], (MJD_to_datetime(float(toa.TOA)) - MJD_to_datetime(bat[i])).seconds, toa.TOAsigma, self.err[i], toalist[i].flags['i'], toalist[i].flags['EMAX']
                print i+1, 'TOA', toalist[i+1].TOA, bat[i+1], (MJD_to_datetime(float(toalist[i+1].TOA)) - MJD_to_datetime(bat[i+1])).seconds, toalist[i+1].TOAsigma, self.err[i+1], toalist[i+1].flags['i'],toalist[i].flags['EMAX']

        for name in self.groups.keys():
            self._wrms(name)



    def tempo2fit(self, toafile=None):
        """
        model.tempo2fit(toafile, pulsefile=None):
        use tempo2 to fit the parfile and toafile.
        """
        self.START = toafile.start
        self.FINISH = toafile.end
        if toafile == None:
            toafile = ''
        if self.TempoVersion == 1:
            self.convert(2)
            self.write()
        if toafile.issubgroup or not toafile.format == 'Tempo2':
            tmptimfile = tempfile.NamedTemporaryFile(mode='w', suffix=".tim", delete=False).name
            tf = open(tmptimfile, 'w')
            tf.write(toafile.tempo2fmt())
            tf.close()
            toafile.toafile = tmptimfile
        self.toafile = toafile
        if toafile.issubgroup:
            jumpdict = {}
            firstgrp = toafile.newjumpgroup[0]
            if firstgrp == 0:
                self.JUMP_0 = (Decimal(0), Decimal(0))
            self.jumps = {}
            for newjp in toafile.newjumpgroup.keys():
                if not newjp == 0:
                    oldjp = toafile.newjumpgroup[newjp]
                    jumpdict['JUMP -jump %d' % newjp] = (self.__dict__['JUMP -jump %d' % oldjp][0] - self.__dict__['JUMP -jump %d' % firstgrp][0], self.__dict__['JUMP -jump %d' % oldjp][1] + self.__dict__['JUMP -jump %d' % firstgrp][1] )
                    self.jumps[newjp] = jumpdict['JUMP -jump %d' % newjp][0]
            if firstgrp == 0:
                del self.__dict__['JUMP -jump 0'] 
            for key in self.__dict__.keys():
                if not key.find('JUMP') == -1:
                    del self.__dict__[key]
                    self.manifest.pop(self.manifest.index(key))
                    del self.parameters[key] #clean up old jumps
            for key in jumpdict.keys():
                self.__dict__[key] = jumpdict[key]
                self.manifest.append(key)
                self.parameters[key] = '1' #default to be fitted

        #take care of DMX ranges outside and inside the [START, FINISH]
        if self.__dict__.has_key('DMX') or self.__dict__.has_key('DMX_0001'):
            DMX, DMXErr, DMXR1, DMXR2 = self.dmxlist   
            for key in DMXR1:
                if DMXR2[key] < self.START or DMXR1[key] > self.FINISH:
                    self.parameters['DMX_' + str(key).rjust(4,'0')] = '0'
                else:
                    pass
                #self.parameters['DMX_' + str(key).rjust(4,'0')] = '1'

        self.groups = toafile.groups.copy()
        self.jumpgroups = toafile.jumpgroups.copy()
        self.jumpgroups.update({0:Decimal(0)}) #set the 0 group to match the jumps initialized in PARfile
        self.write()
        if not len(self.jumps) == len(self.jumpgroups) - 1:
            print "Number of jumps in parfile (%s) doesn't match the number of jumps in the toafile (%s)" % (len(self.jumps), len(self.jumpgroups) - 1)
            #print "This could be a problem."
        if not set(self.jumps.keys()) < set(self.jumpgroups.keys()):
            print "Jump group flag mismatch! try othergroups"
            if set(self.jumps.keys()) - set(toafile.othergroups.keys()) == set([0]):
                print 'othergroups matches'
                for key in set(self.jumps.keys())&set(toafile.othergroups.keys()):
                    self.jumpgroups[key] = toafile.othergroups[key]
                    toafile.jumpgroups[key] = toafile.othergroups[key]
            else:
                print set(self.jumps.keys()) - set(toafile.othergroups.keys())
                print "Groups defined in parfile and not in toafile: %s" % (set(self.jumps.keys()) - set(self.jumpgroups.keys()))
                print "Groups defined in toafile and not in parfile: %s" % (set(self.jumpgroups.keys()) - set(self.jumps.keys()))

        line = getoutput("tempo2 -f %s %s -tempo1 -output general2 -s '{bat} {freq} {post} {err} {post_phase} {binphase} {pre} {pre_phase} {npulse}\n' -outfile resid.dat -showchisq > tmp2log1" % (self.parfile, toafile.toafile))
        line = getoutput("tempo2 -f %s %s -tempo1 -showchisq -newpar > tmp2log2; cat tmp2log2 | grep 'Chisqr/nfree'" % (self.parfile, toafile.toafile))
        self.line = line
        chisqlist = []
        for l in line.split('\n'):
            l  = l.split()[6]
            chisq = l.split('/')[0]
            dof = int(l.split('/')[1])
            try:
                chisq = float(chisq)
            except:
                chisq = 9999999.
            chisqlist.append(chisq)
        #print chisq, dof
        self.chisq = min(chisqlist)
        self.dof = dof
        self.newpar = PARfile('new.par')

        from numpy import genfromtxt
        data = genfromtxt('resid.dat')
        self.toa = np.array(data[...,0])
        self.freq = np.array(data[...,1])
        self.res = np.array(data[...,2]) * 1.e6
        self.err = np.array(data[...,3]) 
        self.phase = np.array(data[...,4])
        self.ophase = np.array(data[...,5])
        self.prefitres = np.array(data[...,6]) * 1.e6
        self.prefitphase = np.array(data[...,7])
        self.npulse = np.array(data[...,8])
        self.weight = 1./self.err**2
        #self.toa = readcol('resid.dat', 0)
        #self.freq = readcol('resid.dat', 1)
        #self.res = readcol('resid.dat', 2)
        #self.err = readcol('resid.dat', 3)*1.e-6
        #self.postphase = readcol('resid.dat', 4)
        from datatools.MJD import MJD_to_datetime
        from datetime import timedelta
        bat = self.toa
        toalist = self.toafile.toalist
        for i in range(len(bat)):
            toa = toalist[i]
            if MJD_to_datetime(float(toa.TOA)) - MJD_to_datetime(bat[i])  > timedelta(seconds = 600):
                for m in range(30,1,-1):
                    print i-m, 'TOA', toalist[i-m].TOA, bat[i-m], (MJD_to_datetime(float(toalist[i-m].TOA)) - MJD_to_datetime(bat[i-m])).seconds, toalist[i-m].TOAsigma, self.err[i-m], toalist[i-m].flags['i']
                print '>', i, 'TOA', toa.TOA, bat[i], (MJD_to_datetime(float(toa.TOA)) - MJD_to_datetime(bat[i])).seconds, toa.TOAsigma, self.err[i], toalist[i].flags['i']
                print i+1, 'TOA', toalist[i+1].TOA, bat[i+1], (MJD_to_datetime(float(toalist[i+1].TOA)) - MJD_to_datetime(bat[i+1])).seconds, toalist[i+1].TOAsigma, self.err[i+1], toalist[i+1].flags['i']

        for name in self.groups.keys():
            self._wrms(name)


    def _wrms(self, *xargs): #must be run after a fit
        from numpy import mean, sum, nan
        if len(xargs) == 0:
            idx = range(len(self.res))
        else:
            idx = self.groups[xargs[0]] #take a particular group
        res = self.res[idx]
        err = self.err[idx]
        weight = 1./err**2
        wmres = sum(res*weight)/sum(weight)
        mres = mean(res)
        sumres = sum(res)
        wsum = sum(weight)
        try:
            wrms = sqrt(sum(res**2*weight - wmres*sumres)/wsum)
            #wrms = sqrt(sum(res**2*weight/wsum - wmres**2))
        except ValueError:
            #print sum(res**2*weight - wmres*sumres)
            #print wsum, wmres, sumres
            #print res
            #print ValueError
            wrms = nan
        if len(xargs) == 0:
            self.wrms = {'all':wrms}
        else:
            name = xargs[0] 
            if not self.__dict__.has_key('wrms'):
                self._wrms()
            self.wrms[name] = wrms


    def plotres(self, *xargs):
        """The same as the tempo/util/plotres, use as model.plotres(), can take general matplotlib plotting parameters as keywords."""
        import datatools.MJD
        import pylab as plt
        zero = np.array([0]*len(self.toa))
        if len(xargs) == 0:
            TOAdates = [datatools.MJD.MJD_to_datetime(t) for t in self.toa]
            plt.plot(TOAdates, zero, '--')
            plt.errorbar(TOAdates, self.res, yerr = self.err, fmt = '.')
            plt.xlabel('Date')
            plt.ylabel(r'residual (${\rm \mu}$s)')
        else:
            if xargs[0] == 'freq':
                plt.plot(self.freq, zero, '--')
                plt.errorbar(self.freq, self.res, yerr = self.err, fmt = '.')
                plt.xlim(600., 3100.)
                plt.xlabel('Frequency (MHz)')
                plt.ylabel(r'residual (${\rm \mu}$s)')
            elif xargs[0] == 'ophase':
                plt.plot(self.ophase, zero, '--')
                plt.errorbar(self.ophase, self.res, yerr = self.err, fmt = '.')
                plt.xlim(0., 1.)
                plt.xlabel('Orbital Phase')
                plt.ylabel(r'residual (${\rm \mu}$s)')
            else:pass
        plt.show()

    def plot(self, Xlabel, Ylabel, groups=None, colors=None, ax=None, fig=None, LegendOn=False, LegendLabels=None, **kwargs):
        """ploting routine, 
        possible X-axis: number, date, mjd, ophase, freq 
        possilbe Y-axis: err, res, averes, prefit, DMX 
            """
        c = colors
        from pylab import subplot, xlabel, ylabel, legend
        from datatools.MJD import MJD_to_datetime
        from numpy import inf, nan
        colors = c
        labeldict = {
                "number":"Number",
                "date":"Date",
                "year":"day of the year",
                "mjd":"MJD",
                "phase":'phase',
                "err":'error',
                "error":'error',
                "freq":"Frequency (MHz)",
                "freqency":"Frequency (MHz)",
                "ophase":"Orbital Phase",
                "res":r'residual (${\rm \mu}$s)',
                "residual":r'residual (${\rm \mu}$s)',
                "averes":r'averaged residual (${\rm \mu}$s)',
                "aveerr":r'averaged erro (${\rm \mu}$s)',
                "post phase":'postfit phase',
                "prefit":r'prefit residual (${\rm \mu}$s)',
                "DMX":r'Delta DM (pc cm$^{-3}$)',
                }
        def dayofyear(dt):
            return float(dt.strftime('%j')) + float(dt.strftime('%H'))/24
        if ax == None:
            ax = subplot(111)
        if Ylabel == "DMX":
            from matplotlib.ticker import FormatStrFormatter
            majorFormatter = FormatStrFormatter('%5.0f')
            DMX, DMXErr, DMXR1, DMXR2 = self.dmxlist
            DMXR = []
            DMXRErr = []
            DMXvalue = []
            DMXerror = []
            for i in DMX.keys():
                if Xlabel == "date":
                    DMXR.append(MJD_to_datetime(float(DMXR1[i]+DMXR2[i])/2))
                    DMXRErr.append((MJD_to_datetime(float(DMXR2[i])) - MJD_to_datetime(float(DMXR1[i])))/2)
                elif Xlabel == "year":
                    DMXR.append(dayofyear(MJD_to_datetime(float(DMXR1[i]+DMXR2[i])/2)))
                    DMXRErr.append((MJD_to_datetime(float(DMXR2[i])) - MJD_to_datetime(float(DMXR1[i])))/2)
                else:
                    DMXR.append((float(DMXR1[i]+DMXR2[i])/2))
                    DMXRErr.append((float(DMXR2[i]) - float(DMXR1[i]))/2)
                    ax.xaxis.set_major_formatter(majorFormatter)
                DMXvalue.append(float(str(DMX[i])))
                DMXerror.append(float(str(DMXErr[i])))
            try:
                if colors == None:
                    ax.errorbar(DMXR, DMXvalue, xerr=DMXRErr, yerr=DMXerror, fmt='.', **kwargs)
                else:
                    ax.errorbar(DMXR, DMXvalue, xerr=DMXRErr, yerr=DMXe, fmt='.', color=colors[grp], **kwargs)
            except ValueError:
                print 'R: ', DMXR
                print 'value: ', DMXvalue
                print 'Rerr: ', DMXRErr
                print 'DM err: ', DMXerror
            xlabel(labeldict[Xlabel])
            ylabel(labeldict[Ylabel])
            return ax
        if groups == None and not Ylabel in ['averes','aveerr']:
            groups = self.groups.keys()
            groups.sort()
            #groups = self.toafile.grouporder
        elif groups == None and Ylabel in ['averes', 'aveerr']:
            groups = self.averes.keys()
            groups.sort()
        if Xlabel == "date":
            Xlimit = [MJD_to_datetime(100000), MJD_to_datetime(0.)]
        else:
            Xlimit = [inf,0]
        subplots = {}
        for grp in groups:
            Yerr = None
            if Ylabel == "averes":
                if Xlabel == "date":
                    X = [MJD_to_datetime(t) for t in self.avetoa[grp]]
                elif Xlabel == 'year':
                    X = [dayofyear(MJD_to_datetime(t)) for t in self.avetoa[grp]]
                elif Xlabel == 'mjd':
                    from matplotlib.ticker import FormatStrFormatter
                    majorFormatter = FormatStrFormatter('%5.0f')
                    ax.xaxis.set_major_formatter(majorFormatter)
                    X = [float(t) for t in self.avetoa[grp]]
                elif Xlabel == "ophase":
                    raise "Not allowed to plot averes on phase"
                elif Xlabel == "freq" or Xlabel == "frequency":
                    raise "Not allowed to plot averes vs frequence"
                elif Xlabel == "number":
                    raise "Not allowed to plot averes vs number"
                elif Xlabel == "phase":
                    raise "Not allowed to plot averes vs phase"
                else:
                    raise "X label %s not allowed" % Xlabel
                Y = self.averes[grp]
                Yerr = self.aveerr[grp]
            else:
                idx = self.groups[grp]
                if Ylabel == "res" or Ylabel == "residual":
                    Y = self.res[idx]
                    Yerr = self.err[idx]
                elif Ylabel == "err" or Ylabel == "error":
                    Y = self.err[idx]
                elif Ylabel == "prefit":
                    Y = self.prefitres[idx]
                    Yerr = self.err[idx]
                elif Ylabel == "post phase":
                    Y = self.phase[idx]
                if Xlabel == "date":
                    X = [MJD_to_datetime(t) for t in self.toa[idx]]
                elif Xlabel == 'year':
                    X = [dayofyear(MJD_to_datetime(t)) for t in self.toa[idx]]
                    ax.set_xlim([0,366])
                elif Xlabel == "mjd":
                    from matplotlib.ticker import FormatStrFormatter
                    majorFormatter = FormatStrFormatter('%5.0f')
                    ax.xaxis.set_major_formatter(majorFormatter)
                    X = [float(t) for t in self.toa[idx]]
                elif Xlabel == "ophase":
                    X = self.ophase[idx]
                elif Xlabel == "freq" or Xlabel == "frequency":
                    X = self.freq[idx]
                elif Xlabel == "number":
                    X = idx
                elif Xlabel == "phase":
                    X = self.phase[idx]
                else:
                    raise "X label %s not allowed" % Xlabel
            if colors == None:
                if Yerr == None:
                    subp = ax.plot(X, Y, '.', markeredgewidth=0, label=grp, **kwargs)
                else:
                    subp = ax.errorbar(X, Y, Yerr, fmt='.', mew=0, label=grp, **kwargs)
            else:
                if Yerr == None:
                    subp = ax.plot(X, Y, '.', markeredgewidth=0, color=colors[grp], label=grp, **kwargs)
                else:
                    #subp = ax.errorbar(X, Y, Yerr, fmt='.', color=colors[grp], markeredgewidth=1, mec=colors[grp], label=grp, **kwargs)
                    subp = ax.errorbar(X, Y, Yerr, fmt='.', color=colors[grp], mec=colors[grp], label=grp, **kwargs)
            subplots.update({grp:subp})
            Xlimit[0] = min(min(X), Xlimit[0])
            Xlimit[1] = max(max(X), Xlimit[1])
        #ax.plot(Xlimit,[0.,0.],'k--')
        ax.axhline(y=0, color='k', linestyle='--')
        ax.set_xlabel(labeldict[Xlabel])
        ax.set_ylabel(labeldict[Ylabel])
        if LegendOn:
            if LegendLabels == None:
                legend([subplots[g] for g in groups], [g for g in groups], loc=2, numpoints=1)
            else:
                legend([subplots[g] for g in groups], LegendLabels, loc=2, numpoints=1)
        return ax                


    def average(self, groups='', lapse=0.5):
        """
        Performing the daily-average, 
        arguments: groups = [list of groups one wants to average] (Default all groups)
                   lapse=0.5 (the window size for daily averaging, default to half day)
        """
        from numpy import mean, sum, array
        toafile = self.toafile
        toa = self.toa
        res = self.res
        err = self.err
        weight = self.weight
        self.avetoa = {}
        self.averes = {}
        self.aveerr = {}
        self.avewrms = {}
        self.mediansigma = {}
        self.toagrps = {}
        self.toagrpkeys = {}
        if groups == '':
            keys = self.groups.keys()
        elif groups == 'allinone':
            self.groups.update({'all':range(len(toafile.toalist))})
            keys = ['all']
        else:
            keys = groups
        for key in keys:
            idx = self.groups[key]
            toagrp = [(toa[i], i) for i in idx]
            toagrp.sort(key = lambda x: x[0])
            grpkey = toagrp[0][0]
            subgrp = {}
            subgrp[grpkey] = [toagrp[0][1]]
            grpkeylist =[grpkey]
            self.avetoa[key] = []
            self.averes[key] = []
            self.aveerr[key] = []
            for i in [x[1] for x in toagrp[1:]]:
                t = toa[i]
                if abs(t - grpkey) <= lapse:
                    subgrp[grpkey].append(i)
                else:
                    grpkey = t
                    grpkeylist.append(grpkey)
                    subgrp[grpkey] = [i]
            self.toagrps[key]  = subgrp
            self.toagrpkeys[key] = grpkeylist
            for grpkey in grpkeylist:
                idx = subgrp[grpkey]
                avetoa = mean(toa[idx])
                weightsum = sum(weight[idx])
                ressum = sum(res[idx]*weight[idx])
                averes = ressum/weightsum
                N = len(idx)
                if N > 1:
                    sigma = sqrt((sum([(res[i] - averes)**2*weight[i] for i in idx])/(N-1))/weightsum)
                else:
                    sigma = err[idx[0]]
                if sigma == 0.:
                    print "zero res err caused by similar residual in the same group -- ignored"
                    print idx
                    print res[idx], averes
                    print toa[idx]
                    for t in [ self.toafile.toalist[i] for i in idx]:
                        print t.file, t.TOA
                        print t.line
                self.avetoa[key].append(avetoa) 
                self.averes[key].append(averes)
                self.aveerr[key].append(sigma)
            self.avetoa[key] = array(self.avetoa[key])
            self.averes[key] = array(self.averes[key])
            self.aveerr[key] = array(self.aveerr[key])
            def _wrms(res, err):
                nres = []
                nerr = []
                for i in range(len(err)):
                    if not err[i] == 0.:
                        nres.append(res[i])
                        nerr.append(err[i])
                res = np.array(nres)
                err = np.array(nerr)
                weight = 1./err**2
                wmres = sum(res*weight)/sum(weight)
                mres = np.mean(res)
                sumres = sum(res)
                wsum = sum(weight)
                try:
                    wrms = sqrt(sum(res**2*weight - wmres*sumres)/wsum)
                except ValueError:
                    return None
                return wrms
            self.avewrms[key] = _wrms(self.averes[key],self.aveerr[key])
            self.mediansigma[key] = np.median(self.aveerr[key])

def parseerror(val, err):
    if not type(val) in [Decimal, __Decimal]:
        val = __Decimal(str(val))
    if not type(err) in [Decimal, __Decimal]:
        err = __Decimar(str(err))
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


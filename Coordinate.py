import re
from numpy import float64 as float
from decimal import Decimal

def SSprec(SS, error=9.):
    p = re.compile('\d+\.*(?P<digits>\d*)')
    m = p.search(str(SS))
    #print m.group('digits')
    return error*0.1**len(m.group('digits'))

class RA(object):
    def __init__(self, value):
        if isinstance(value, (tuple, list)):
            if value[1] == 'None':
                error = 9
            else:
                error = value[1]
            value = value[0]
            try:
                self.in_unit_degree = float(value)
                self.uncertainty = SSprec(value, error=float(error))
                hh = self.in_unit_degree/15
                self.HH = int(hh)
                mm = (hh % 1)*60
                self.MM = int(mm)
                self.SS = (mm % 1)*60
                self.SSuncertainty = self.uncertainty*60.*4.
            except ValueError:
                pattern = re.compile('(?P<HH>\d\d)[:\s](?P<MM>\d\d)[:\s](?P<SS>\d\d\.*\d*)')
                match = pattern.search(value)
                self.HH = int(match.group('HH'))
                self.MM = int(match.group('MM'))
                self.SS = float(match.group('SS'))
                self.SSuncertainty = SSprec(match.group('SS'), error=float(error))
                self.in_unit_degree = (self.HH+self.MM/60.+self.SS/60./60.)*15
                self.uncertainty = self.SSuncertainty/60./4.
        elif isinstance(value, (basestring,float,Decimal)):
            return RA.__init__(self, (value, 9))
    def __str__(self):
        return '%s:%s:%s' % (self.HH,self.MM,self.SS)
class Dec(object):
    def __init__(self, value):
        if isinstance(value, (tuple, list)):
            if value[1] == 'None':
                error = 9
            else:
                error = value[1]
            value = value[0]
            try:
                self.in_unit_degree = float(value)
                self.uncertainty = SSprec(value, error=float(error))
                if self.in_unit_degree < 0.:
                    self.sign = '-'
                else:
                    self.sign = '+'
                dd = abs(self.in_unit_degree)
                self.dd = int(dd)
                mm = (dd % 1)*60
                self.mm = int(mm)
                self.ss = (mm % 1)*60
                self.SSuncertainty = self.uncertainty*60.*60
            except ValueError:
                pattern = re.compile('(?P<sign>-*\+*)(?P<dd>\d\d)[:\s](?P<mm>\d\d)[:\s](?P<ss>\d\d\.*\d*)')
                match = pattern.search(value)
                self.sign = str(match.group('sign'))
                self.dd = int(match.group('dd'))
                self.mm = int(match.group('mm'))
                self.ss = float(match.group('ss'))
                self.SSuncertainty = SSprec(match.group('ss'), error=float(error))
                if self.sign == '-':
                    self.in_unit_degree = -1.*(self.dd+self.mm/60.+self.ss/60./60.)
                else:
                    self.in_unit_degree = self.dd+self.mm/60.+self.ss/60./60.
                self.uncertainty = self.SSuncertainty/60./60.
                #print self.SSuncertainty
        elif isinstance(value, (basestring,float,Decimal)):
            return Dec.__init__(self, (value, 9))
    def __str__(self):
        return '%s%s:%s:%s' % (self.sign,self.dd,self.mm,self.ss)



class coordinate(object):
    def __init__(self, values):
        if isinstance(values, basestring):
            ra,dec = values.split('  ')[:2]
            self.RA = RA(ra)
            self.Dec = Dec(dec)
        elif isinstance(values, tuple) and len(values) == 2:
            self.RA = RA(str(values[0]))
            self.Dec = Dec(str(values[1]))
        else:
            raise 'Unrecognized Input'
    def __str__(self):
        return str(str(self.RA)+ '  ' + str(self.Dec))

class ecs(object):
    '''Equatorial coordinate system'''
    def __init__(self, RA, Dec):
        self.RA = RA
        self.Dec = Dec
    def __str__(self):
        return '%s(%g S) %s(%g s)' % (self.RA, self.RA.SSuncertainty ,self.Dec,self.Dec.SSuncertainty)
    def ecs2gcs(self):
        ra = self.RA.in_unit_degree
        dec = self.Dec.in_unit_degree
        GCOra = 266.405100
        GCOdec = -28.936175
        GCNPra = 192.859508
        GCNPdec = 27.128336


class gcs(object):
    '''Galactic coordinate system'''
    def __init__(self, l, b):
        self.l = l
        self.b = b
    def __str__(self):
        return '%s(%g S) %s(%g s)' % (self.l, self.l.SSuncertainty ,self.b,self.b.SSuncertainty)

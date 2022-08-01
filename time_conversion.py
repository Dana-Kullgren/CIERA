

from numpy import *
import gzip
import math
import subprocess

def read_units(filestring):
    """reads from the supplied conv file and stores the physical units"""
    s = filestring+'.conv.sh'
    f=open(s,'r')
    count=0
    while count<10:
        a=f.readline()
        b=a.split('=')
        if b[0]=='massunitcgs':
            massunitcgs = float(b[1])
            count+=1
        elif b[0]=='massunitmsun':
            massunitmsun = float(b[1])
            count+=1
        elif b[0]=='mstarunitcgs':
            mstarunitcgs = float(b[1])
            count+=1
        elif b[0]=='mstarunitmsun':
            mstarunitmsun = float(b[1])
            count+=1
        elif b[0]=='lengthunitcgs':
            lengthunitcgs = float(b[1])
            count+=1
        elif b[0]=='lengthunitparsec':
            lengthunitparsec = float(b[1])
            count+=1
        elif b[0]=='timeunitcgs':
            timeunitcgs = float(b[1])
            count+=1
        elif b[0]=='timeunitsmyr':
            timeunitsmyr = float(b[1])
            count+=1
        elif b[0]=='nbtimeunitcgs':
            nbtimeunitcgs = float(b[1])
            count+=1
        elif b[0]=='nbtimeunitsmyr':
            nbtimeunitsmyr = float(b[1])
            count+=1
    f.close()
    units = []
    unittype = [('m_cgs', float), ('m_msun', float), ('mstar_cgs', float), ('mstar_msun', float), ('l_cgs', float), ('l_pc', float), ('t_cgs', float),('t_myr', float), ('nbt_cgs', float), ('nbt_myr', float)]
    units.append((massunitcgs, massunitmsun, mstarunitcgs, mstarunitmsun, lengthunitcgs, lengthunitparsec, timeunitcgs, timeunitsmyr, nbtimeunitcgs, nbtimeunitsmyr))
    units = array(units, dtype=unittype)
    return units






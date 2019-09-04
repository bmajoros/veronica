#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Author: William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import ProgramName
from Rex import Rex
rex=Rex()
from Interval import Interval

def loadRepeats(filename):
    repeats=set()
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=2): continue
            repeats.add(int(fields[0]))
    return repeats

def parseInterval(s):
    if(not rex.find("\((\d+),(\d+)\)",s)): 
        raise Exception("can't parse interval: "+s)
    return Interval(int(rex[1]),int(rex[2]))

def loadProcessed(filename):
    junctions=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split("\t")
            if(len(fields)!=7): continue
            if(fields[4]!="EXON_DELETED"): continue
            leftInterval=parseInterval(fields[5])
            rightInterval=parseInterval(fields[6])
            junctions.append([leftInterval,rightInterval])
    return junctions

def distanceToRepeat(interval,repeats):
    smallest=None
    for x in repeats:
        d=interval.distanceFromPoint(x)
        if(smallest is None or d<smallest): smallest=d
    return smallest
        
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <repeat-positions.txt> <postprocessed.txt>\n")
(repeatFile,processedFile)=sys.argv[1:]

# Load data
repeats=loadRepeats(repeatFile)
junctions=loadProcessed(processedFile)

# Process junctions
for junction in junctions:
    (leftInterval,rightInterval)=junction
    D1=distanceToRepeat(leftInterval,repeats)
    D2=distanceToRepeat(rightInterval,repeats)
    print(D1,D2,sep="\t")







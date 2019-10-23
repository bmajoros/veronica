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

def loadTargetSites(filename):
    bySite={}
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=3): continue
            (site,pos,seq)=fields
            bySite[site]=int(pos)
    return bySite

def getDistance(site1,site2,bySite):
    pos1=bySite[site1]
    pos2=bySite[site2]
    dist=abs(pos2-pos1)
    return dist

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <postprocessed.txt> <target-sites.txt>\n")
(postprocessedFile,targetsFile)=sys.argv[1:]

sites=loadTargetSites(targetsFile)
with open(postprocessedFile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=6): continue
        (PAIR,normCount,factor,site1,site2,rawCount)=fields
        dist=getDistance(site1,site2,sites)
        print(dist,normCount,sep="\t")



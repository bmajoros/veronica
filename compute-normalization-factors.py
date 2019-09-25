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

class Site:
    def __init__(self,ID,pos):
        self.ID=ID
        self.pos=pos
        self.norm=None

def loadSites(filename):
    sites=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=3): continue
            (siteID,pos,seq)=fields
            pos=int(pos)
            site=Site(siteID,pos)
            sites.append(site)
    sites.sort(key=lambda x: x.pos)
    return sites

def processDepth(filename,sites):
    # Precondition: the input file and sites list are both sorted by position
    siteIndex=0
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=3): continue
            (ID,pos,depth)=fields
            pos=int(pos); depth=int(depth)
            while(siteIndex<len(sites) and pos>sites[siteIndex].pos):
                siteIndex+=1
            if(siteIndex>=len(sites)): break
            site=sites[siteIndex]
            if(pos==site.pos): site.norm=depth

def maxDepth(sites):
    m=None
    for site in sites:
        if(m is None or site.norm>m): m=site.norm
    return m

def normalize(sites):
    m=maxDepth(sites)
    for site in sites:
        site.norm/=m

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <nolib.depth> <target-sites.txt>\n")
(depthFile,sitesFile)=sys.argv[1:]

# Load data
sites=loadSites(sitesFile)
processDepth(depthFile,sites)

# Normalize counts into a pseudo-probability
normalize(sites)

# Generate output
for site in sites:
    print(site.ID,site.pos,site.norm,sep="\t")


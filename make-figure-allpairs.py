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

BASE="/data/gersbachlab/bill"
TARGET_SITES=BASE+"target-sites.txt"

class Record:
    def __init__(self,read,leftTuple,rightTuple,strand,deleted):
        self.read=read
        self.leftTuple=leftTuple
        self.rightTuple=rightTuple
        self.strand=strand
        self.deleted=deleted
    def getGuideDistLen(self):
        (leftGuide,leftDist,leftLen)=parseTuple(self.leftTuple)
        (rightGuide,rightDist,rightLen)=parseTuple(self.rightTuple)
        return (leftGuide,leftDist,leftLen,rightGuide,rightDist,rightLen)

def getExtents(records,targetSites):
    smallest=None; largest=None
    for rec in records:
        (leftGuide,leftDist,leftLen,rightGuide,rightDist,rightLen)=\
            rec.getGuideDistLen()
        leftTarget=targetSites.get(leftGuide,None)
        rightTarget=targetSites.get(rightGuide,None)
        if(leftTarget is None or rightTarget is None):
            raise Exception("Guide not found")
        leftBegin=leftTarget+leftDist
        rightBegin=rightTarget+rightDist

def loadTargetSites(filename):
    h={}
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=3): continue
            (guide,pos,seq)=fields
            h[guide]=int(pos)
    return h

def parseDist(dist):
    if(not rex.find("D=([\d-]+)",dist)): raise Exception("Can't parse: "+dist)
    return int(rex[1])

def parseLen(L):
    if(not rex.find("L=([\d-]+)",L)): raise Exception("Can't parse: "+L)
    return int(rex[1])

def parseTuple(tup):
    fields=tup.split()
    if(len(fields)!=3): raise Exception("Can't parse: "+tup)
    (guide,dist,length)=fields
    dist=parseDist(dist)
    length=parseLen(length)
    return (guide,dist,length)

def incCounts(begin,end,counts):
    for i in range(begin,end): counts[i]+=1

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <postprocessed.txt>\n")
(infile,)=sys.argv[1:]

# Load the target site locations
targetSites=loadTargetSites(TARGET_SITES)

# Load the edit data
records=[]
with open(infile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split("\t")
        if(len(fields)!=5): continue
        (read,leftTuple,rightTuple,strand,deleted)=fields
        if(deleted!="EXON_DELETED"): continue
        rec=Record(read,leftTuple,rightTuple,strand,deleted)
        records.append(rec)
(BEGIN,seqLen)=getExtents(records,targetSites)

# Compile pileups
counts=[0]*seqLen
for rec in records:
        (leftGuide,leftDist,leftLen,rightGuide,rightDist,rightLen)=\
            rec.getGuideDistLen()
        leftTarget=targetSites.get(leftGuide,None)
        rightTarget=targetSites.get(rightGuide,None)
        if(leftTarget is None or rightTarget is None):
            raise Exception("Guide not found")
        leftBegin=leftTarget+leftDist
        rightBegin=rightTarget+rightDist
        adjLeft=leftBegin-BEGIN
        adjRight=rightBegin-BEGIN
        incCounts(adjLeft-leftLen,adjLeft,counts)
        incCounts(adjRight,adjRight+rightLen,counts)

# Print out the pileups
for i in range(len(counts)):
    if(counts[i]>0): print(i,counts[i],sep="\t")


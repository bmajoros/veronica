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

BASE="/data/gersbachlab/bill"
TARGET_SITES=BASE+"/target-sites.txt"

class Record:
    def __init__(self,read,leftTuple,rightTuple,strand,deleted,int1,int2):
        self.read=read
        self.leftTuple=leftTuple
        self.rightTuple=rightTuple
        self.strand=strand
        self.deleted=deleted
        self.interval1=int1
        self.interval2=int2
    def getGuideDistLen(self):
        (leftGuide,leftDist,leftLen)=parseTuple(self.leftTuple)
        (rightGuide,rightDist,rightLen)=parseTuple(self.rightTuple)
        return (leftGuide,leftDist,leftLen,rightGuide,rightDist,rightLen)

def getExtents(records,targetSites):
    smallest=None; largest=None
    for rec in records:
        m=min(rec.interval1.begin,rec.interval1.end,
              rec.interval2.begin,rec.interval2.end)
        if(smallest is None or m<smallest): smallest=m
        m=max(rec.interval1.begin,rec.interval1.end,
              rec.interval2.begin,rec.interval2.end)
        if(largest is None or m>largest): largest=m
    return (smallest,largest)

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

def pickPoint(interval,target):
    return interval.begin if \
        abs(interval.begin-target)<abs(interval.end-target) \
        else interval.end

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
        if(len(fields)!=7): continue
        (read,leftTuple,rightTuple,strand,deleted,interval1,interval2)=fields
        if(deleted!="EXON_DELETED"): continue
        interval1=Interval.parseInt(interval1)
        interval2=Interval.parseInt(interval2)
        rec=Record(read,leftTuple,rightTuple,strand,deleted,interval1,interval2)
        records.append(rec)
(BEGIN,seqLen)=getExtents(records,targetSites)
BEGIN=0

# Compile pileups
counts=[0]*seqLen
for rec in records:
    (leftGuide,leftDist,leftLen)=parseTuple(rec.leftTuple)
    (rightGuide,rightDist,rightLen)=parseTuple(rec.rightTuple)
    leftTarget=targetSites.get(leftGuide,None)
    rightTarget=targetSites.get(rightGuide,None)
    if(leftTarget is None or rightTarget is None):
        raise Exception("Guide not found")
    #incCounts(rec.interval1.begin-BEGIN,rec.interval1.end-BEGIN,counts)
    #incCounts(rec.interval2.begin-BEGIN,rec.interval2.end-BEGIN,counts)
    leftPoint=pickPoint(rec.interval1,leftTarget)
    rightPoint=pickPoint(rec.interval2,rightTarget)
    incCounts(leftPoint,leftPoint+1,counts)
    incCounts(rightPoint,rightPoint+1,counts)

# Print out the pileups
for i in range(len(counts)):
    if(counts[i]>0): print(i,counts[i],sep="\t")


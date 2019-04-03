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
TARGETS=BASE+"target-sites.txt"

LEFT_TARGET=4800 #V_50_9  4800    CACCACTCACCTC
RIGHT_TARGET=6411 #V_51_1  6411    GACCATTTCCCAC
LEFT_GUIDE="V_50_9"
RIGHT_GUIDE="V_51_1"
BEGIN=0 # 3800
END=12000 #7411
SEQ_LEN=END-BEGIN

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

def incCounts_OLD(begin,end,counts):
    for i in range(0,begin): counts[i]+=1
    for i in range(end,SEQ_LEN): counts[i]+=1

def incCounts(begin,end,counts):
    #print(begin,end)
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

counts=[0]*SEQ_LEN
with open(infile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split("\t")
        if(len(fields)!=7): continue
        (read,leftTuple,rightTuple,strand,deleted,interval1,interval2)=fields
        if(deleted!="EXON_DELETED"): continue
        interval1=Interval.parseInt(interval1)
        interval2=Interval.parseInt(interval2)
        (leftGuide,leftDist,leftLen)=parseTuple(leftTuple)
        (rightGuide,rightDist,rightLen)=parseTuple(rightTuple)
        #incCounts(interval1.begin-BEGIN,interval1.end-BEGIN,counts)
        #incCounts(interval2.begin-BEGIN,interval2.end-BEGIN,counts)
        leftPoint=pickPoint(interval1,LEFT_TARGET)
        rightPoint=pickPoint(interval2,RIGHT_TARGET)
        incCounts(leftPoint,leftPoint+1,counts)
        incCounts(rightPoint,rightPoint+1,counts)

for i in range(len(counts)):
    if(counts[i]>0): print(i,counts[i],sep="\t")
    #print(i,0,sep="\t")


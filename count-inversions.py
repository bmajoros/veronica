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

MAX_DIST=25
BASE="/data/gersbachlab/bill"
TARGET_SITES=BASE+"/target-sites.txt"

class Record:
    def __init__(self,read,leftTuple,rightTuple,strand,deleted,int1,int2):
        self.read=read
        self.leftTuple=leftTuple
        self.rightTuple=rightTuple
        (guide,D,L)=self.parseTuple(leftTuple)
        self.leftGuide=guide
        self.leftDist=D
        self.leftLen=L
        (guide,D,L)=self.parseTuple(rightTuple)
        self.rightGuide=guide
        self.rightDist=D
        self.rightLen=L
        self.strand=strand
        self.deleted=deleted
        self.interval1=int1
        self.interval2=int2
    def print(self):
        print(self.read,self.leftTuple,self.rightTuple,self.strand,
              self.deleted,self.interval1.toString(),self.interval2.toString(),
              sep="\t")
    def parseDist(self,dist):
        if(not rex.find("D=([\d-]+)",dist)):
            raise Exception("Can't parse: "+dist)
        return int(rex[1])
    def parseLen(self,L):
        if(not rex.find("L=([\d-]+)",L)): raise Exception("Can't parse: "+L)
        return int(rex[1])
    def parseTuple(self,tup):
        fields=tup.split()
        if(len(fields)!=3): raise Exception("Can't parse: "+tup)
        (guide,dist,length)=fields
        dist=self.parseDist(dist)
        length=self.parseLen(length)
        return (guide,dist,length)

def loadRecords(filename):
    records=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split("\t")
            if(len(fields)!=7): continue
            (read,leftTuple,rightTuple,strand,deleted,interval1,interval2)=\
                fields
            if(deleted!="EXON_DELETED"): continue
            interval1=Interval.parseInt(interval1)
            interval2=Interval.parseInt(interval2)
            rec=Record(read,leftTuple,rightTuple,strand,deleted,interval1,\
                           interval2)
            records.append(rec)
    return records

def loadTargetSites(filename):
    h={}
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=3): continue
            (guide,pos,seq)=fields
            h[guide]=int(pos)
    return h

def findTarget(guide,targetSites):
    target=targetSites.get(guide,None)
    if(target is None): raise Exception("Cannot find: "+guide)
    return target

def intervalDistance(interval,target):
    if(interval.contains(target)): return 0
    d1=abs(target-interval.begin)
    d2=abs(target-interval.end)
    d=min(d1,d2)
    return d

def getSide(interval,target):
    if(interval.begin>=target): return "RIGHT"
    if(interval.end<=target): return "LEFT"
    return None

def classify(rec,targetSites):
    leftTarget=findTarget(rec.leftGuide,targetSites)
    rightTarget=findTarget(rec.rightGuide,targetSites)
    leftSide=getSide(rec.interval1,leftTarget)
    rightSide=getSide(rec.interval2,rightTarget)
    d=intervalDistance(rec.interval1,leftTarget)
    #if(d>MAX_DIST):
    #    rec.print()
    #    raise Exception("Distance is too great: "+str(d)+" "+\
    #                        rec.interval1.toString()+" "+str(leftTarget))
    #d=intervalDistance(rec.interval2,rightTarget)
    #if(d>MAX_DIST):
    #    rec.print()
    #    raise Exception("Distance is too great: "+str(d)+" "+\
    #                        rec.interval2.toString()+" "+str(rightTarget))
    if(leftSide is None or rightSide is None): return None
    if(leftSide=="LEFT" and rightSide=="RIGHT"): return "EXON_DELETION"
    if(leftSide=="LEFT" and rightSide=="LEFT"): return "INVERSION"
    if(leftSide=="RIGHT" and rightSide=="RIGHT"): return "INVERSION"
    if(leftSide=="RIGHT" and rightSide=="LEFT"): return "CIRCULARIZATION"

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <postprocessed.txt>\n")
(filename,)=sys.argv[1:]

# Load the target site locations and alignment results
targetSites=loadTargetSites(TARGET_SITES)
records=loadRecords(filename)

# Process records
for rec in records:
    c=classify(rec,targetSites)
    if(c is None): continue
    print(c)



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
import gzip

MIN_DELETION=100

class Tuple:
    def __init__(self,guide,distance,length):
        self.guide=guide
        self.distance=distance
        self.length=length
    def toString(self):
        return self.guide+" [D="+str(self.distance)+"] L="+\
            str(self.length)

class Record:
    def __init__(self,readID,tuple1,tuple2,deleted,strand,int1,int2):
        self.readID=readID
        self.match1=self.parseTuple(tuple1)
        self.match2=self.parseTuple(tuple2)
        self.deleted=deleted=="EXON_DELETED"
        self.sameStrands=strand=="+"
        self.interval1=self.parseInterval(int1)
        self.interval2=self.parseInterval(int2)
    def getDeletionLen(self):
        L=None
        if(self.interval1.begin<self.interval2.begin):
            L=self.interval2.begin-self.interval1.end
        else:
            L=self.interval1.begin-self.interval2.end
        #print(L)
        return L if L>=0 else 0
    def parseInterval(self,interval):
        if(not rex.find("(\d+):(\d+)",interval)):
            raise Exception("Can't parse interval: "+interval)
        return Interval(int(rex[1]),int(rex[2]))
    def getKey(self):
        return self.match1.toString()+" "+self.match2.toString()
    def print(self):
        match1=self.match1
        match2=self.match2
        if(match1.guide>match2.guide): (match1,match2)=(match2,match1)
        print(self.readID,match1.toString(),match2.toString(),
              "+" if self.sameStrands else "-",
              "EXON_DELETED" if self.deleted else "",
              self.interval1.toString(),self.interval2.toString(),
              sep="\t")
    def parseTuple(self,tuple1):
        if(not rex.find("(\S+)\s+\[D=([\d-]+)\]\s+L=(\d+)",tuple1)):
            raise Exception("Can't parse: "+tuple1)
        return Tuple(rex[1],int(rex[2]),int(rex[3]))

def findBestGuides(records):
    counts={}
    for rec in records:
        if(not rec.deleted): continue
        counts[rec.match1.guide]=counts.get(rec.match1.guide,0)+1
        counts[rec.match2.guide]=counts.get(rec.match2.guide,0)+1
    array=[]
    for key in counts.keys():
        count=counts[key]
        array.append([key,count])
    array.sort(key=lambda x: x[1],reverse=True)
    for elem in array:
        (guide,count)=elem
        print(guide,count,sep="\t")

def findBestExamples(records):
    records.sort(key=
                 lambda x: 1/(1+abs(x.match1.distance)) * 1/(1+abs(x.match2.distance)) * abs(x.match1.length) * abs(x.match2.length) * (1 if x.deleted else 0),
                 reverse=True)
    n=100000
    if(n>len(records)): n=len(records)
    for i in range(n):
        rec=records[i]
        if(rec.getDeletionLen()<MIN_DELETION): continue
        rec.print()

def countDeleted(records):
    n=0; sameStrand=0
    for rec in records:
        if(rec.deleted): 
            n+=1
            if(rec.sameStrands): sameStrand+=1
    return (n,sameStrand)

def countDiffStrand(records):
    n=0
    for rec in records:
        if(not rec.sameStrands): n+=1
    return n

def computeNorm(pairKey,normFactors):
    if(not rex.find("(\S+) (\S+)",pairKey)):
        raise Exception("Can't parse pair key: "+pairKey)
    (id1,id2)=(rex[1],rex[2])
    norm1=normFactors[id1]
    norm2=normFactors[id2]
    norm=1.0/(1-(1-norm1)*(1-norm2))
    return norm

def countGuidePairs(records,normFactors):
    pairCounts={}
    for rec in records:
        if(not rec.deleted): continue
        key=""
        if(rec.match1.guide<rec.match2.guide):
            key=rec.match1.guide+" "+rec.match2.guide
        else:
            key=rec.match2.guide+" "+rec.match1.guide
        pairCounts[key]=pairCounts.get(key,0)+1
    array=[]
    for key in pairCounts.keys():
        count=pairCounts[key]
        norm=computeNorm(key,normFactors)
        count*=norm
        array.append([count,key])
    array.sort(key=lambda x: -x[0])
    for rec in array:
        (count,key)=rec
        norm=computeNorm(key,normFactors)
        print(round(count,1),round(norm,2),key,sep="\t")

def readNormFactors(filename):
    factors={}
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=3): continue
            (ID,pos,P)=fields
            P=float(P)
            factors[ID]=P
    return factors

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <processed.txt> <norm-factors.txt>\n")
(filename,normFile)=sys.argv[1:]

# Read normalization factors
normFactors=readNormFactors(normFile)

# Read records
records=[]
seen=set()
IN=gzip.open(filename,"rt") if rex.find(".gz$",filename) \
    else open(filename,"rt")
#with open(filename,"rt") as IN:
for line in IN:
    fields=line.rstrip("\n").split("\t")
    if(len(fields)!=7): raise Exception(line)
    (readID,tuple1,tuple2,strand,deleted,interval1,interval2)=fields
    rec=Record(readID,tuple1,tuple2,deleted,strand,interval1,interval2)
    key=rec.getKey()
    if(key not in seen): records.append(rec)
    seen.add(key)
IN.close()

# Compute statistics
n=len(records)
(numDel,numDelSameStrand)=countDeleted(records)
percentDel=round(float(numDel)/float(n),3)
percentDelMinus=round(float(numDel-numDelSameStrand)/float(numDel),3)
print("percent deleted: ",percentDel*100,"% = ",numDel,"/",n,", ",
      percentDelMinus*100,"% of which are -",sep="")

diffStrand=countDiffStrand(records)
percentDiffStrand=round(float(diffStrand)/float(n),3)
print("percent different strand: ",percentDiffStrand*100,"% = ",
      diffStrand,"/",n,sep="")

countGuidePairs(records,normFactors)

findBestExamples(records)

findBestGuides(records)

# M03884:303:000000000-C4RM6:1:1111:23549:4101	V_51_28 [D=9] L=90	V_50_29 [D=11] L=60	EXON_DELETED



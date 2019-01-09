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

class Tuple:
    def __init__(self,guide,distance,length):
        self.guide=guide
        self.distance=distance
        self.length=length
    def toString(self):
        return self.guide+" [D="+str(self.distance)+"] L="+\
            str(self.length)

class Record:
    def __init__(self,readID,tuple1,tuple2,deleted,strand):
        self.readID=readID
        self.match1=self.parseTuple(tuple1)
        self.match2=self.parseTuple(tuple2)
        self.deleted=deleted=="EXON_DELETED"
        self.sameStrands=strand=="+"
    def getKey(self):
        return self.match1.toString()+" "+self.match2.toString()
    def print(self):
        match1=self.match1
        match2=self.match2
        if(match1.guide>match2.guide): (match1,match2)=(match2,match1)
        print(self.readID,match1.toString(),match2.toString(),
              "+" if self.sameStrands else "-",
              "EXON_DELETED" if self.deleted else "",sep="\t")

    def parseTuple(self,tuple1):
        if(not rex.find("(\S+)\s+\[D=(\d+)\]\s+L=(\d+)",tuple1)):
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
                 lambda x: 1/(1+x.match1.distance) * 1/(1+x.match2.distance) * \
                     x.match1.length * x.match2.length * \
                     (1 if x.deleted else 0),
                 reverse=True)
    for i in range(100):
        records[i].print()

def countDeleted(records):
    n=0
    for rec in records:
        if(rec.deleted): n+=1
    return n

def countDiffStrand(records):
    n=0
    for rec in records:
        if(not rec.sameStrands): n+=1
    return n

def countGuidePairs(records):
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
        array.append([count,key])
    array.sort(key=lambda x: -x[0])
    for rec in array:
        print(rec[0],rec[1],sep="\t")


#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <outputs.txt>\n")
(filename,)=sys.argv[1:]

# Read records
records=[]
seen=set()
with open(filename,"rt") as IN:
    for line in IN:
        fields=line.rstrip("\n").split("\t")
        if(len(fields)!=5): raise Exception(line)
        (readID,tuple1,tuple2,strand,deleted)=fields
        rec=Record(readID,tuple1,tuple2,deleted,strand)
        key=rec.getKey()
        if(key not in seen): records.append(rec)
        seen.add(key)

# Compute statistics
n=len(records)
numDel=countDeleted(records)
percentDel=round(float(numDel)/float(n),3)
print("percent deleted: ",percentDel*100,"% = ",numDel,"/",n,sep="")

diffStrand=countDiffStrand(records)
percentDiffStrand=round(float(diffStrand)/float(n),3)
print("percent different strand: ",percentDiffStrand*100,"% = ",
      diffStrand,"/",n,sep="")

countGuidePairs(records)

findBestExamples(records)

findBestGuides(records)

# M03884:303:000000000-C4RM6:1:1111:23549:4101	V_51_28 [D=9] L=90	V_50_29 [D=11] L=60	EXON_DELETED



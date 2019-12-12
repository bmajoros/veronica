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

class GuidePair:
    def __init__(self,key,rank,norm,abundance,rawCount):
        self.key=key
        self.rank=rank
        self.normalized=norm
        self.abundance=abundance
        self.rawCount=rawCount

def load(filename):
    recs={}
    with open(filename,"rt") as IN:
        rank=1
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=7): continue
            if(fields[0]!="PAIR"): continue
            (pair,normalizedScore,normFactor,abundance,guide1,guide2,rawCount)=\
                fields
            key=guide1+" "+guide2
            rec=GuidePair(key,rank,normalizedScore,abundance,rawCount)
            recs[key]=rec
            rank+=1
    return recs

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <renormalized1-3.txt> <renormalized4.txt>\n")
(file1,file2)=sys.argv[1:]

hash1=load(file1)
hash2=load(file2)
combined=[]
for key in hash1.keys():
    rec1=hash1[key]
    rec2=hash2.get(key,None)
    if(rec2 is None): 
        #print("NO RECORD FOUND IN "+file2+" for "+key)
        continue
    aveRank=(rec1.rank+rec2.rank)/2
    aveRank=int(aveRank) if aveRank==int(aveRank) else aveRank
    combined.append([aveRank,key,rec1.rawCount,rec2.rawCount,rec1.normalized,
                     rec2.normalized,rec1.abundance,rec1.rank,rec2.rank])
combined.sort(key=lambda x: x[0])
for rec in combined:
    (aveRank,key,rawCount1,rawCount2,normalized1,normalized2,abundance,
     rank1,rank2)=rec
    print(aveRank,rank1,rank2,key,rawCount1,rawCount2,normalized1,normalized2,
          abundance,sep="\t")


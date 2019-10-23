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

STRAND_FILE="/home/bmajoros/charlie/veronica/long-guide-sequences.txt"

def loadStrand(filename):
    byGuide={}
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=7): continue
            (guide,chrom,pos1,pos2,strand,intron,seq)=fields
            byGuide[guide]=strand
    return byGuide

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <postprocessed.txt>\n")
(filename,)=sys.argv[1:]

strands=loadStrand(STRAND_FILE)
with open(filename,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        #if(len(fields)!=5 or fields[0]!="GUIDE"): continue
        #(GUIDE,normCount,factor,guide,raw)=fields
        #strand=strands[guide]
        #strand=1 if strand=="+" else -1
        #print(strand,normCount,sep="\t")
        if(len(fields)!=6 or fields[0]!="PAIR"): continue
        (PAIR,normCount,factor,guide1,guide2,raw)=fields
        strand1=strands[guide1]
        strand2=strands[guide2]
        strand1=1 if strand1=="+" else -1
        strand2=1 if strand2=="+" else -1
        strand=strand1*strand2
        print(strand,normCount,sep="\t")




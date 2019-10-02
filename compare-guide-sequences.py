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

GUIDES_FILE="/home/bmajoros/charlie/veronica/long-guide-sequences.txt"

def loadGuides(filename):
    guides={}
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)<7): continue
            (ID,chrom,begin,end,strand,intron,seq)=fields
            #print(ID,seq.upper(),sep="\t")
            guides[seq]=ID
    return guides

def countMismatches(str1,str2):
    n=0
    for x,y in zip(str1,str2):
        if(x!=y): n+=1
    return n

#=========================================================================
# main()
#=========================================================================
guides=loadGuides(GUIDES_FILE)
guideStrings=[x for x in guides.keys()]
numGuides=len(guideStrings)
for i in range(numGuides):
    guideI=guideStrings[i]
    for j in range(numGuides):
        if(i==j): continue
        guideJ=guideStrings[j]
        mis=countMismatches(guideI,guideJ)
        print(mis,guides[guideI],guides[guideJ],sep="\t")



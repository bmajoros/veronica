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
from Guide import Guide
from GuidePair import GuidePair
from Deletion import Deletion

def parsePair(fields,pairs):
    #  fields:  1   V_50_51   V_51_39
    pair=GuidePair(fields[0],fields[1],fields[2])
    if(pair.guide1.intron!=50 or pair.guide2.intron!=51):
        raise Exception("Unexpected exons: "+guide1+" && "+guide2)
    pairs.append(pair)

def parseDeletion(fields,deletions):
    # fields:  K00282:378:H3T5TBBXY:1:1113:19806:45203 V_50_26 [D=0] L=76  V_51_38 [D=0] L=75   +   EXON_DELETED  (3142,3218)  (11083,11158)
    (readID,guide1,d1,l1,guide2,d2,l2,sign,deleted,coords1,coords2)=fields
    d=Deletion(readID,guide1,d1,l1,guide2,d2,l2,sign,deleted,coords1,coords2)
    deletions.append(d)

def parseGuideCount(fields,guides):
    # V_50_18 48
    guide=Guide(fields[0])
    count=int(fields[1])
    guides.append([guide,count])

def parseFile(filename,pairs,guides,deletions):
    with open(filename,"rt") as IN:
        IN.readline(); IN.readline()
        for line in IN:
            fields=line.rstrip().split()
            numFields=len(fields)
            if(numFields==3): parsePair(fields,pairs)
            elif(numFields==2): parseGuideCount(fields,guides)
            elif(numFields==11): parseDeletion(fields,deletions)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <postprocessed.txt>\n")
(infile,)=sys.argv[1:]

pairs=[] # array of Pair objects
guideCounts=[] # array of pairs: (Guide,int)
deletions=[] # array of Deletion objects
parseFile(infile,pairs,guideCounts,deletions)

# How many guides?
numGuides=len(guideCounts)
print(numGuides,"guides")
closeGuides=set(); exactGuides=set()
for d in deletions:
    if(abs(d.d1)<=10): closeGuides.add(d.guide1.toString())
    if(abs(d.d2)<=10): closeGuides.add(d.guide2.toString())
    if(abs(d.d1)==0): exactGuides.add(d.guide1.toString())
    if(abs(d.d2)==0): exactGuides.add(d.guide2.toString())
numCloseGuides=len(closeGuides)
print(numCloseGuides,"guides cut within 10 bp")
numExactGuides=len(exactGuides)
print(numExactGuides,"guides cut exactly")
for guideCount in guideCounts:
    (guide,count)=guideCount
    if(guide.toString() not in exactGuides):
        print("\t"+guide.toString(),"was never cut exactly")

# How many guide pairs?
numPairs=len(pairs)
print(numPairs,"pairs found")

# How many pairs had a perfect cut at both target sites?
perfectDeletions=set()
for d in deletions:
    if(d.d1==0 and d.d2==0):
        key=d.guide1.toString()+" "+d.guide2.toString()
        perfectDeletions.add(key)
numPerfect=len(perfectDeletions)
print(numPerfect,"guide pairs with perfect deletions")

# How many pairs had a nearly perfect cut at both target sites?
closeDeletions=set()
for d in deletions:
    if(abs(d.d1)<=10 and abs(d.d2)<=10):
        key=d.guide1.toString()+" "+d.guide2.toString()
        closeDeletions.add(key)
        #print("CLOSE:\t"+d.toString())
numClose=len(closeDeletions)
print(numClose,"guide pairs with nearly perfect deletions (within 10 bp)")


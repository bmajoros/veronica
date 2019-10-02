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
from FastqReader import FastqReader
from Rex import Rex
rex=Rex()
from Translation import Translation

MAX_MISMATCH=5 # No two guides differ by fewer than 7 positions
GUIDES_FILE="/home/bmajoros/charlie/veronica/long-guide-sequences.txt"
GUIDE1_BEGIN=30
GUIDE2_BEGIN=125

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

def fuzzyMatch(guideSeq,guideStrings,guides):
    bestGuide=None; bestMismatch=None
    for guide in guideStrings:
        mis=countMismatches(guideSeq,guide)
        if(bestMismatch is None or mis<bestMismatch):
            bestGuide=guide
            bestMismatch=mis
    if(bestMismatch>MAX_MISMATCH): bestGuide=None
    return None if bestGuide is None else guides[bestGuide]

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <infile.fastq.gz> <min-quality:30>\n")
(infile,MIN_QUALITY)=sys.argv[1:]
MIN_QUALITY=int(MIN_QUALITY)

shouldRev=rex.find("R2",infile)
begin=GUIDE1_BEGIN if not shouldRev else GUIDE2_BEGIN
end=begin+21
guides=loadGuides(GUIDES_FILE)
guideStrings=guides.keys()
reader=FastqReader(infile)
while(True):
    rec=reader.nextSequence()
    if(rec is None): break
    (readID,seq,qual,rawQual,pair)=rec
    guideSeq=seq[begin:end]
    qualSeq=qual[begin:end]
    rawQual=rawQual[begin:end]
    if(shouldRev): guideSeq=Translation.reverseComplement(guideSeq)
    if(min(qualSeq)<MIN_QUALITY): continue
    if(len(guideSeq)!=21): raise Exception("Guide is wrong length")
    ID=guides.get(guideSeq,None)
    if(ID is None): ID=fuzzyMatch(guideSeq,guideStrings,guides)
    #if(ID is not None): print(ID,readID,sep="\t")
    score=round(sum(qualSeq)/len(qualSeq),2)
    if(ID is not None): print(ID,score,sep="\t")


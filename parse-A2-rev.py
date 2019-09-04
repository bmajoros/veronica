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
from Pipe import Pipe
from Shuffler import Shuffler

SHUFFLE=False
MIN_MATCH=14 # 1% of permuted sequences have 14 matches or more
BASE="/home/bmajoros/charlie/veronica/newdata"
FWD_READS=BASE+"/fastq/A2_S1_L001_R1_001.fastq.gz"
REV_READS=BASE+"/fastq/A2_S1_L001_R2_001.fastq.gz"
GUIDES=BASE+"/guides.txt"
ALIGNER="/home/bmajoros/BOOM/ungapped-aligner"
MATRIX="/home/bmajoros/alignment/matrices/NUC.4.4"

#R1: tatatatcttgtggaaaggacgaaacaccg (then 21 bp intron 50 gRNA) gttttagtactctggaaacagaatctactaaaacaaggcaaaatgccgtgtttatctcgtcaacttgttggcgagatttttttGTACACCGGTGCCCA

# Chop off 30 bp, then take the next 21 bp as the guide; use ungapped-aligner to match it to a guide, allowing mismatches

#R2: cactaaagggaacaaaagctggagctccaccgcgggcggccgcaaaaaaatctcgccaacaagttgacgagataaacacggcattttgccttgttttagtagattctgtttccagagtactaaaac (21 bp intron 51 gRNA) ca

# Chop off 126 bp, then take 21 bp guide, use ungapped-aligner to match to a guide

def loadGuides(filename):
    guides50=[]; guides51=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=2): continue
            (ID,seq)=fields
            seq=seq.upper()
            if(SHUFFLE): seq=Shuffler.shuffleString(seq)
            guide=[ID,seq]
            if(not rex.find("V_(\d+)_\d+",ID)):
                raise Exception("Can't parse guide ID"+ID)
            if(rex[1]=="50"): guides50.append(guide)
            else: guides51.append(guide)
    return (guides50,guides51)

def revGuides(guides):
    revs=[]
    for guide in guides:
        (ID,seq)=guide
        rev=Translation.reverseComplement(seq)
        revs.append([ID,rev])
    return revs

def matchGuide(guideSeq,guides):
    #for guide in guides:
    #    (ID,seq)=guide
    #    if(guideSeq==seq): return (guide,21)
    bestGuides=None; bestScore=None
    for guide in guides:
        (ID,seq)=guide
        score=align(guideSeq,seq)
        if(bestScore is None or score>bestScore):
            bestGuide=guide; bestScore=score
    if(bestScore<MIN_MATCH): return None
    return (bestGuide,bestScore)

def align(a,b):
    L=len(a)
    matches=0
    for i in range(L):
        if(a[i]==b[i]): matches+=1
    return matches

#=========================================================================
# main()
#=========================================================================
(guides50,guides51)=loadGuides(GUIDES)
rev50=revGuides(guides50); rev51=revGuides(guides51)
#guides50.extend(rev50); guides51.extend(rev51)
reader=FastqReader(REV_READS)
#numFound=0
#n=0
begin=126
while(True):
    #n+=1
    #if(n>=100000): break ### TESTING
    rec=reader.nextSequence()
    if(rec is None): break
    (ID,seq,qual,pair)=rec
    guideSeq=seq[begin:(begin+21)]
    match=matchGuide(guideSeq,guides51)
    #if(match is None): match=matchGuide(guideSeq,rev50)
    if(match is None): continue
    (bestMatch,score)=match
    #print(guideSeq,bestMatch,score,sep="\t")
    (bestID,bestSeq)=bestMatch
    print(ID,bestID,score,sep="\t")
    #numFound+=1
reader.close()
#print(numFound,"found")


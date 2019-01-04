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
from SamReader import SamReader
from FastaReader import FastaReader

MIN_SOFT_MASK=16
MAX_DISTANCE=25
GENOME="/home/bmajoros/veronica/exon51.fasta"

class Target:
    def __init__(self,ID,pos,seq):
        self.ID=ID
        self.pos=pos
        self.seq=seq

def softMaskLength(cigar):
    for op in cigar.ops:
        if(op.op=="S"): return op.length
    return None

def goodCigar(cigar):
    L=cigar.length()
    if(L!=2): return False
    s=softMaskLength(cigar)
    if(s==None or s<MIN_SOFT_MASK): return False
    return True

def findBreakpoint(rec):
    return rec.refPos+rec.CIGAR[0].length

def loadTargets(filename):
    targets=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=3): continue
            (ID,pos,seq)=fields
            targets.append(Target(ID,int(pos),seq))
    return targets

def findTarget(targets,breakpoint):
    bestTarget=None
    bestDiff=None
    for target in targets:
        diff=abs(breakpoint-target.pos)
        if(bestDiff is None or diff<bestDiff):
            bestTarget=target
            bestDiff=diff
    return bestTarget

def getUnaligned(rec):
    seq=rec.seq
    cigar=rec.CIGAR
    if(cigar[0].op=="S"):
        return seq[:cigar[0].length]
    if(cigar[1].op=="S"):
        return seq[cigar[0].length:]
    else: raise Exception("internal error")

def findUnaligned(unaligned,genome):
    L=len(genome)
    querylen=len(unaligned)
    last=L-querylen
    for i in range(last):
        if(genome[i:(i+querylen)]==unaligned):
            return i

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <filename.sam> <target-sites.txt>\n")
(samFile,targetFile)=sys.argv[1:]

# Load genomic sequence
(Def,genome)=FastaReader.firstSequence(GENOME)

# Load target locations
targets=loadTargets(targetFile)

# Process SAM file
reader=SamReader(samFile)
while(True):
    rec=reader.nextSequence()
    if(rec is None): break
    if(rec.flag_unmapped()): continue
    if(rec.CIGAR.completeMatch()): continue
    if(not goodCigar(rec.CIGAR)): continue
    #print(rec.ID,rec.CIGAR.toString(),sep="\t")
    breakpoint=findBreakpoint(rec)
    nearestTarget=findTarget(targets,breakpoint)
    distance=abs(breakpoint-nearestTarget.pos)
    if(distance>MAX_DISTANCE): continue

    # Try to align the unaligned part to the other intron
    unaligned=getUnaligned(rec)
    pos=findUnaligned(unaligned,genome)
    print(breakpoint,pos,sep="\t")
    #print(rec.ID,unaligned,sep="\t")

    #print(nearestTarget.pos,breakpoint,sep="\t")
    #print(rec.ID,rec.CIGAR.toString(),
          #"mult="+str(rec.flag_hasMultipleSegments()),
          #"aligned="+str(rec.flag_properlyAligned()),
          #"unmapped="+str(rec.flag_unmapped()),
          #"first="+str(rec.flag_firstOfPair()),
          #"second="+str(rec.flag_secondOfPair()),
          #"rev="+str(rec.flag_revComp()),
          #"secondary="+str(rec.flag_secondaryAlignment()),
          #"failed="+str(rec.flag_failedFilters()),
          #"dup="+str(rec.flag_PCRduplicate()),
          #"suppl="+str(rec.flag_supplAlignment()),
          #sep="\t")




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
import os
import ProgramName
from SamReader import SamReader
from FastaReader import FastaReader
from FastaWriter import FastaWriter
from Translation import Translation
from Interval import Interval
from Pipe import Pipe
from Rex import Rex
rex=Rex()
import TempFilename
from CigarString import CigarString

BOOM="/home/bmajoros/BOOM"
MATRIX="/home/bmajoros/alignment/matrices/NUC.4.4"
GAP_OPEN=2
GAP_EXTEND=1
MIN_SOFT_MASK=16
MIN_MATCH=16
MAX_DISTANCE=25
GENOME="/home/bmajoros/veronica/exon51.fasta"
TARGETS="/home/bmajoros/veronica/target-sites.txt"
FASTA_WRITER=FastaWriter()

class Target:
    def __init__(self,ID,pos,seq):
        self.ID=ID
        self.pos=pos
        self.seq=seq
        if(not rex.find("V_(\d+)_\d+",ID)):
            raise Exception("Can't parse target ID")
        self.intron=rex[1]

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
    matchLen=0
    for op in rec.CIGAR.ops:
        if(op.getOp()=="M"): matchLen=op.getLength()
    begin=rec.refPos+rec.CIGAR[0].length
    return (begin,matchLen)

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
        #print("XXX",seq[:cigar[0].length],sep="\t")
        return seq[:cigar[0].length]
    if(cigar[0].op=="M" and cigar[1].op=="S"):
        #print("YYY",seq[cigar[0].length:],sep="\t")
        return seq[cigar[0].length:]
    else: raise Exception("internal error")

def longestMatchFromCigar(cigar):
    L=cigar.length()
    longestOpIndex=None
    longestOpLen=0
    for i in range(L):
        op=cigar[i]
        if(op.getOp()!="M"): continue
        opLen=op.getLength()
        if(opLen>longestOpLen):
            longestOpLen=opLen
            longestOpIndex=i
    if(longestOpIndex is None): return None
    readPos=0; genomePos=0
    for i in range(longestOpIndex):
        op=cigar[i]
        opLen=op.getLength()
        if(op.advanceInQuery()): readPos+=opLen
        if(op.advanceInRef()): genomePos+=opLen
    return [readPos,genomePos,longestOpLen]

def findUnaligned(unaligned,genome,revGenome):
    #print("Running SmithWaterman on:",unaligned)
    cigar=CigarString(smithWaterman(unaligned,genome,GAP_OPEN,GAP_EXTEND))
    print("SmithWaterman: ",cigar.toString())
    longest=longestMatchFromCigar(cigar)
    if(longest is None): longest=(-1,-1,-1)
    (readPos,genomePos,matchLen)=longest
    strand="+"
    if(matchLen<MIN_MATCH):
        unaligned=Translation.reverseComplement(unaligned)
        cigar=CigarString(smithWaterman(unaligned,genome,GAP_OPEN,
                                         GAP_EXTEND))
        longest=longestMatchFromCigar(cigar)
        if(longest is None): return None
        (readPos,genomePos,matchLen)=longest
        readPos=len(unaligned)-readPos-1 ###
        if(matchLen<MIN_MATCH): return None
        strand="-"
    if(readPos<0): raise Exception("error")
    readSubseq=unaligned[readPos:(readPos+matchLen)]
    genomeSubseq=genome[genomePos:(genomePos+matchLen)]
    return (genomePos,strand,matchLen,genomeSubseq)

def getAlignedIntervals(rec):
    readPos=0
    genomePos=rec.refPos
    cigar=rec.CIGAR
    for i in range(cigar.length()):
        op=cigar[i]
        if(op.op=="M"):
            readInterval=Interval(readPos,readPos+op.length)
            genomeInterval=Interval(genomePos,genomePos+op.length)
            return [readInterval,genomeInterval]
        if(op.advanceInQuery()): readPos+=op.length
        if(op.advanceInRef()): genomePos+=op.length
    return None

def sanityCheckAlignment(rec,genome):
    intervals=getAlignedIntervals(rec)
    (readPos,genomePos)=intervals
    readSeq=rec.getSequence()[readPos.getBegin():readPos.getEnd()]
    genomeSeq=genome[genomePos.getBegin():genomePos.getEnd()]
    return (readPos,genomePos,readSeq,genomeSeq)

def getMatchProportion(seq1,seq2):
    nMatch=0
    nMismatch=0
    for i in range(len(seq1)):
        if(seq1[i]==seq2[i]): nMatch+=1
        else: nMismatch+=1
    return float(nMatch)/float(nMatch+nMismatch)

def writeFile(defline,seq):
    filename=TempFilename.generate("fasta")
    FASTA_WRITER.writeFasta(defline,seq,filename)
    return filename

def swapInsDel(cigar):
    # This is done because my aligner defines insertions and deletions
    # opposite to how they're defined in the SAM specification
    newCigar=""
    for x in cigar:
        if(x=="I"): x="D"
        elif(x=="D"): x="I"
        newCigar+=x
    return newCigar

def smithWaterman(seq1,seq2,gapOpen,gapExtend):
    file1=writeFile("read",seq1)
    file2=writeFile("genome",seq2)
    cmd=BOOM+"/smith-waterman -q "+MATRIX+" "+str(gapOpen)+" "+\
        str(gapExtend)+" "+file1+" "+file2+" DNA"
    output=Pipe.run(cmd)
    #print(cmd); print(output); #exit()
    os.remove(file1)
    os.remove(file2)
    if(not rex.find("CIGAR=(\S+)",output)):
        raise Exception("Can't parse aligner output: "+output)
    cigar=rex[1]
    cigar=swapInsDel(cigar) # because my aligner defines cigars differently
    return cigar

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <read-ID> <sam-file>\n")
(readID,samFile)=sys.argv[1:]

targets=loadTargets(TARGETS)

# Load genomic sequence
(Def,genome)=FastaReader.firstSequence(GENOME)
revGenome=Translation.reverseComplement(genome)
reader=SamReader(samFile)
while(True):
    rec=reader.nextSequence()
    if(rec is None): break
    if(rec.ID!=readID): continue
    if(rec.flag_unmapped()): continue
    if(rec.CIGAR.completeMatch()): continue
    if(not goodCigar(rec.CIGAR)): continue
    print("RAW READ =",rec.seq)
    (breakpoint,anchorLen1)=findBreakpoint(rec)
    nearestTarget1=findTarget(targets,breakpoint)
    distance1=abs(breakpoint-nearestTarget1.pos)
    if(distance1>MAX_DISTANCE): continue
    print("Bowtie:",rec.CIGAR.toString())
    (readPos,genomePos,readSeq,genomeSeq1)=sanityCheckAlignment(rec,genome)
    print("G1=",genomePos,"\tD1=",distance1,"\tL1=",anchorLen1,"\t",\
              genomeSeq1,sep="")

    # Try to align the unaligned part to the other intron
    unaligned=getUnaligned(rec)
    unalignedPos=findUnaligned(unaligned,genome,revGenome)
    if(unalignedPos is None): continue
    (pos,strand2,anchorLen2,genomeSeq2)=unalignedPos
    nearestTarget2=findTarget(targets,pos)
    distance2=abs(pos-nearestTarget2.pos)
    if(distance2>MAX_DISTANCE): continue
    exonDeleted="EXON_DELETED" if nearestTarget1.intron!=nearestTarget2.intron\
        else ""
    print("G2=(",pos,",",pos+anchorLen2,")\tD2=",distance2,"\tL2=",anchorLen2,\
              "\t",strand2,"\t",genomeSeq2,sep="")
    #print(rec.getID(),"\t",
    #      nearestTarget1.ID," [D=",distance1,"] L=",anchorLen1,"\t",
    #      nearestTarget2.ID," [D=",distance2,"] L=",anchorLen2,"\t",
    #      strand2,"\t",exonDeleted,
    #      sep="")



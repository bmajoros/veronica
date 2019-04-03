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

CHECK_ALIGNMENTS=False
BOOM="/home/bmajoros/BOOM"
MATRIX="/home/bmajoros/alignment/matrices/NUC.4.4"
GAP_OPEN=2
GAP_EXTEND=1
MIN_SOFT_MASK=16 #16
MIN_MATCH=16 #16 is good, 10 is bad
MAX_DISTANCE=25
GENOME="/home/bmajoros/charlie/veronica/exon51.fasta"
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

def longestSoftMask(cigar):
    L=cigar.length()
    longestOpIndex=None
    longestOpLen=0
    for i in range(L):
        op=cigar[i]
        if(op.getOp()!="S"): continue
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

def findUnaligned(unaligned,genome):
    cigar=CigarString(smithWaterman(unaligned,genome,GAP_OPEN,GAP_EXTEND))
    longest=longestMatchFromCigar(cigar)
    if(longest is None): longest=(-1,-1,-1)
    (readPos,genomePos,matchLen)=longest
    strand="+"
    readSubseq=None; genomeSubseq=None
    if(matchLen<MIN_MATCH):
        unaligned=Translation.reverseComplement(unaligned)
        cigar=CigarString(smithWaterman(unaligned,genome,GAP_OPEN,
                                         GAP_EXTEND))
        longest=longestMatchFromCigar(cigar)
        if(longest is None): return None
        (readPos,genomePos,matchLen)=longest
        if(CHECK_ALIGNMENTS):
            readSubseq=unaligned[readPos:(readPos+matchLen)]
            genomeSubseq=genome[genomePos:(genomePos+matchLen)]
        readPos=len(unaligned)-readPos-1 ###
        if(matchLen<MIN_MATCH): return None
        strand="-"
    else:
        if(CHECK_ALIGNMENTS):
            readSubseq=unaligned[readPos:(readPos+matchLen)]
            genomeSubseq=genome[genomePos:(genomePos+matchLen)]
    if(readPos<0): raise Exception("error")
    if(CHECK_ALIGNMENTS):
        matchProp=getMatchProportion(readSubseq,genomeSubseq)
        if(matchProp<0.8):
            print("MISMATCH IN UNALIGNED PORTION:",strand,matchProp)
            print(cigar.toString())
            print(readSubseq+"\n"+genomeSubseq+"\n==========================")
            exit()
    return (genomePos,strand,matchLen,readPos,genomeSubseq,readSubseq)

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
    if(readSeq!=genomeSeq and 
       getMatchProportion(readSeq,genomeSeq)<0.9):
        print("MAIN READ ALIGNMENT")
        print(readSeq+"\n"+genomeSeq+"\n=============================")
        raise Exception("bad alignment")
        return False
    return True

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
    os.remove(file1)
    os.remove(file2)
    if(not rex.find("CIGAR=(\S+)",output)):
        raise Exception("Can't parse aligner output: "+output)
    cigar=rex[1]
    cigar=swapInsDel(cigar) # because my aligner defines cigars differently
    return cigar

def getBowtieMatch(rec):
    cigar=rec.getCigar()
    L=cigar.length()
    for i in range(L):
        if(cigar[i].getOp()=="M"):
            return cigar[i]

def getAllSoftmasks(rec,minLen):
    cigar=rec.getCigar()
    L=cigar.length()
    softmasks=[]
    for i in range(L):
        if(cigar[i].getOp()!="S"): continue
        if(cigar[i].getLength()<minLen): continue
        softmasks.append(cigar[i])
    return softmasks

def getMatchBreakpoint(match,softMask):
    if(softMask.interval1.getEnd()<=match.interval1.getBegin()):
        return (match.interval1.getBegin(),match.interval2.getBegin())
    else:
        return (match.interval1.getEnd(),match.interval2.getEnd())

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
readsSeen=set()
reader=SamReader(samFile)
while(True):
    rec=reader.nextSequence()
    if(rec is None): break
    if(rec.ID in readsSeen): continue
    if(rec.flag_unmapped()): continue
    if(rec.CIGAR.completeMatch()): continue
    rec.CIGAR.computeIntervals(rec.getRefPos())
    readLen=len(rec.getSequence())
    match1=getBowtieMatch(rec)
    anchorLen1=match1.interval1.getLength()
    if(anchorLen1<MIN_MATCH): continue
    softmasks=getAllSoftmasks(rec,MIN_MATCH)
    if(len(softmasks)<1 or len(softmasks)>2): continue
    bestSoftmask=None; bestDistance1=None; bestDistance2=None
    bestReadPos1=None; bestRefPos1=None
    bestAnchorLen2=None; bestReadPos2=None; bestRefPos2=None
    for softmask in softmasks:
        breakpoint1=getMatchBreakpoint(match1,softmask)
        (breakQuery1,breakRef1)=breakpoint1
        nearestTarget1=findTarget(targets,breakRef1)
        distance1=abs(breakRef1-nearestTarget1.pos)
        if(distance1>MAX_DISTANCE): continue
        unaligned=rec.seq[softmask.interval1.begin:softmask.interval1.end]
        unalignedPos=findUnaligned(unaligned,genome)
        if(unalignedPos is None): continue
        (leftPos,strand2,anchorLen2,leftReadPos,genomeSubseq,readSubseq)=\
            unalignedPos
        if(anchorLen2<MIN_MATCH): continue
        leftReadPos+=softmask.interval1.begin
        rightPos=leftPos+anchorLen2
        rightReadPos=leftReadPos+anchorLen2
        nearestTarget2left=findTarget(targets,leftPos)
        distance2left=leftPos-nearestTarget2left.pos
        nearestTarget2right=findTarget(targets,rightPos)
        distance2right=rightPos-nearestTarget2right.pos
        nearestTarget2=None; distances=None; breakRef2=None; readPos=None
        if(softmask.interval2.begin>=match1.interval2.begin):
            nearestTarget2=nearestTarget2left
            distance2=distance2left
            breakRef2=leftPos
            readPos=leftReadPos
        else:
            nearestTarget2=nearestTarget2right
            distance2=distance2right
            breakRef2=rightPos
            readPos=rightReadPos
        if(distance2>MAX_DISTANCE): continue
        if(bestSoftmask is None or 
           anchorLen2>bestAnchorLen2):
            bestSoftmask=softmask; bestAnchorLen2=anchorLen2
            bestDistance1=distance1; bestDistance2=distance2
            bestReadPos1=breakQuery1; bestRefPos1=breakRef1
            bestReadPos2=readPos; bestRefPos2=breakRef2
    if(bestSoftmask is None): continue
    if(False):
        print("=======================================")
        print(rec.ID)
        print("Readlen =",readLen)
        print("Bowtie match = ",match1.interval1.toString(),\
                  match1.interval2.toString(),sep="")
        print("Best softmask= ",bestSoftmask.interval1.toString(),\
                  bestSoftmask.interval2.toString(),sep="")
        print("BREAKPOINT1   = ",breakQuery1," : ",breakRef1," dist=",
              bestDistance1," len=",anchorLen1," ",nearestTarget1.ID,sep="")
        print("BREAKPOINT2   = ",breakRef2," (",strand2,") "," dist=",
              bestDistance2," len=",anchorLen2," ",nearestTarget2.ID,sep="")
        print("DELETION LENGTH = ",abs(breakRef2-breakRef1))
        sanityCheckAlignment(rec,genome) ###
    concordant="CONCORDANT" if strand2=="+" else "DISCORDANT"
    exonDeleted=""; 
    if(nearestTarget1.intron!=nearestTarget2.intron):
        exonDeleted="EXON_DELETED"
    elif(strand2=="-"):
        exonDeleted="?"
    else:
        refDelta=abs(bestRefPos2-bestRefPos1)
        queryDelta=abs(bestReadPos2-bestReadPos1)
        indelLen=queryDelta-refDelta
        if(indelLen==0): continue
        exonDeleted="INDEL:"+str(indelLen)
    interval1=match1.interval2
    interval2=Interval(leftPos,leftPos+anchorLen2)
    if(interval1.begin>interval2.begin):
        (interval1,interval2)=(interval2,interval1)
    print(rec.getID(),"\t",
          nearestTarget1.ID," [D=",bestDistance1,"] L=",anchorLen1,"\t",
          nearestTarget2.ID," [D=",bestDistance2,"] L=",anchorLen2,"\t",
          strand2,"\t",exonDeleted,"\t",
          interval1.begin,":",interval1.end,"\t",
          interval2.begin,":",interval2.end,
          sep="")
    readsSeen.add(rec.ID)




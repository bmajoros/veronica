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
from FastaReader import FastaReader
from Rex import Rex
rex=Rex()

BASE="/home/bmajoros/veronica"
GUIDES=BASE+"/guides.txt"
GENOME=BASE+"/exon51.fasta"

def loadGuides(filename):
    array=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=2): continue
            (ID,seq)=fields
            seq=seq.upper()
            leftRight=None
            if(rex.find("V_50",ID)): leftRight="right"
            elif(rex.find("V_51",ID)): leftRight="left"
            if(leftRight is None): raise Exception("Can't parse left/right")
            array.append([ID,seq,leftRight])
    return array

#=========================================================================
# main()
#=========================================================================

(defline,genome)=FastaReader.firstSequence(GENOME)
L=len(genome)
guides=loadGuides(GUIDES)
for guide in guides:
    (ID,seq,leftRight)=guide
    glen=len(seq)
    last=L-glen
    for i in range(last):
        if(genome[i:(i+glen)]==seq):
            pos=i if leftRight=="left" else i+glen
            print(ID,pos,seq,sep="\t")


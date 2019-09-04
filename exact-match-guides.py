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
from Translation import Translation
from Rex import Rex
rex=Rex()
import gzip

BASE="/home/bmajoros/charlie/veronica/newdata"
GUIDES=BASE+"/guides.txt"

def loadGuides(filename):
    guides50=[]; guides51=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=2): continue
            (ID,seq)=fields
            seq=seq.upper()
            guide=[ID,seq]
            if(not rex.find("V_(\d+)_\d+",ID)):
                raise Exception("Can't parse guide ID"+ID)
            if(rex[1]=="50"): guides50.append(guide)
            else: guides51.append(guide)
    return (guides50,guides51)

def report(line1,intron,guide):
    fields=line1.split()
    readID=fields[0]
    print(readID,intron,guide,sep="\t")

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <*.fastq>\n")
(fastqFile,)=sys.argv[1:]

(guides50,guides51)=loadGuides(GUIDES)

with gzip.open(fastqFile,"rt") as IN:
    while(True):
        line1=IN.readline()
        if(line1 is None): break
        read=IN.readline(); plus=IN.readline(); qual=IN.readline()
        #print(seq)
        rev=Translation.reverseComplement(read)
        found=False
        for guide in guides50:
            if(guide[1] in read or guide[1] in rev):
                report(line1,"V_50",guide[0])
                found=True
                break
        if(found): continue
        for guide in guides51:
            if(guide[1] in read or guide[1] in rev):
                report(line1,"V_51",guide[0])
                found=True
                break
        #if(not found): print("NOT FOUND")

        



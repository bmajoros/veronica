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

BASE="/home/bmajoros/charlie/veronica/newdata"
GUIDEPAIRS=BASE+"/guides.txt"

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

#=========================================================================
# main()
#=========================================================================
(guides50,guides51)=loadGuides(GUIDEPAIRS)
for guide1 in guides50:
    (ID1,seq1)=guide1
    for guide2 in guides51:
        (ID2,seq2)=guide2
        #seq2=Translation.reverseComplement(seq2)
        concatID=ID1+"-"+ID2
        concatSeq=seq1+seq2
        print(">"+concatID+"\n"+concatSeq)



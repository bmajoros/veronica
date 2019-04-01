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

BASE="/home/bmajoros/charlie/veronica/onepair"
TARGETS=BASE+"/../target-sites.txt"

def loadTargets(filename):
    h={}
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=3): continue
            (id,pos,seq)=fields
            h[id]=pos
    return h

def parseTarget(target):
    fields=target.split()
    if(len(fields)!=3): raise Exception("can't parse: "+target)
    id=fields[0]
    if(not rex.find("D=(\d+)",fields[1])):
        raise Exception("can't parse: "+fields[1])
    D=int(rex[1])
    if(not rex.find("L=(\d+)",fields[2])):
        raise Exception("can't parse: "+fields[2])
    L=int(rex[1])
    return (id,D,L)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <postprocessed.txt>\n")
(infile,)=sys.argv[1:]

targets=loadTargets(TARGETS)
with open(infile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split("\t")
        if(len(fields)!=5): continue
        (read,left,right,strand,deleted)=fields
        if(deleted!="EXON_DELETED"): continue
        (leftID,leftD,leftL)=parseTarget(left)
        (rightID,rightD,rightL)=parseTarget(right)
        leftPos=targets[leftID]
        rightPos=targets[rightID]

# M03884:334:000000000-CC4DY:1:1108:22185:9289    V_50_9 [D=16] L=16      V_51_3 [D=24] L=56      +       EXON_DELETED


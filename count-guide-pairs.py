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

def readFile(filename):
    IDs={}
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=2): continue
            (guideID,readID)=fields
            IDs[readID]=guideID
    return IDs

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <file1.txt> <file2.txt>\n")
(file1,file2)=sys.argv[1:]

IDs=readFile(file1)
pairs={}
with open(file2,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=2): continue
        (id2,readID)=fields
        id1=IDs.get(readID,None)
        if(id1 is not None and id2!=id1):
            key=id1+" "+id2
            pairs[key]=pairs.get(key,0)+1
            #print(id1,id2,sep="\t")
keys=[x for x in pairs.keys()]
keys.sort(key=lambda x: -pairs[x])
for key in keys:
    print(key,pairs[key],sep="\t")



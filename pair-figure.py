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

MIN_COUNT=2

def loadFile(filename):
    points={}
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=2): continue
            (x,y)=(int(fields[0]),int(fields[1]))
            if(y>=MIN_COUNT): points[x]=y
    return points

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <file1> <file2>\n")
(filename1,filename2)=sys.argv[1:]

points1=loadFile(filename1)
points2=loadFile(filename2)
for x in points2.keys():
    if(points1.get(x,0)>0):
        print(x,points1[x],points2[x],sep="\t")






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

def loadFile(filename):
    points=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=2): continue
            (x,y)=(int(fields[0]),int(fields[1]))
            points.append([x,y])
    return points

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <with-guides.txt> <nolib.txt> <mask-width> <min-count>\n")
(filename1,filename2,bandwidth,minCount)=sys.argv[1:]
bandwidth=int(bandwidth)
minCount=int(minCount)

withGuides=loadFile(filename1)
nolib=loadFile(filename2)
inNolib=set()
for point in nolib:
    x=point[0]
    for i in range(x-bandwidth,x+bandwidth+1):
        inNolib.add(i)
for point in withGuides:
    (x,y)=point
    if(y>=minCount and x not in inNolib):
        print(x,y,sep="\t")






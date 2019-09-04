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

BASE="/home/bmajoros/charlie/veronica/newdata"
GUIDEPAIRS=BASE+"/guides.txt"

def loadGuides(filename):
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=2): continue
            (ID,seq)=fields
            seq=seq.upper()
            print(">"+ID+"\n"+seq)

#=========================================================================
# main()
#=========================================================================
loadGuides(GUIDEPAIRS)





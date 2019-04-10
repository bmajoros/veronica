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
from FastaReader import FastaReader

BASE="/home/bmajoros/charlie/veronica"
MASKED=BASE+"/exon51.fasta.masked"

#=========================================================================
# main()
#=========================================================================

(defline,seq)=FastaReader.firstSequence(MASKED)
x=0
for c in seq:
    #code=1 if c=="N" else 0
    #print(x,code,sep="\t")
    if(c=="N"): print(x,1,sep="\t")
    x+=1




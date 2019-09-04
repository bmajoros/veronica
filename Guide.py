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
from Rex import Rex
rex=Rex()

# Represents a guide identifier of the form V_50_1, where 50 is the intron 
# number, and 1 is the identifier of the guide in the intron
class Guide:
    def __init__(self,guide):
        (intron,index)=self.parse(guide)
        self.intron=intron
        self.index=index
    def parse(self,guide):
        if(not rex.find("V_(\d+)_(\d+)",guide)): 
            raise Exception("Can't parse guide: "+guide)
        intron=int(rex[1])
        index=int(rex[2])
        return (intron,index)
    def toString(self):
        return "V_"+str(self.intron)+"_"+str(self.index)


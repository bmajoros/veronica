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
from Guide import Guide

# Represents a pair of guide identifiers of the form V_50_1, and their count
class GuidePair:
    def __init__(self,count,guideID1,guideID2):
        self.count=int(count)
        self.guide1=Guide(guideID1)
        self.guide2=Guide(guideID2)
    def toString(self):
        return self.guide1.toString()+" "+self.guide2.toString+" "+\
            str(self.count)




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
from Guide import Guide
from Interval import Interval

# Represents a line from the postprocessing output:
#   K00282:378:H3T5TBBXY:1:1113:19806:45203 V_50_26 [D=0] L=76      V_51_38 [D=0] L=75      +       EXON_DELETED    (3142,3218)     (11083,11158)
class Deletion:
    def __init__(self,readID,guide1,d1,l1,guide2,d2,l2,sign,deleted,
                 coords1,coords2):
        self.readID=readID
        self.guide1=Guide(guide1)
        self.d1=self.parseD(d1)
        self.l1=self.parseL(l1)
        self.guide2=Guide(guide2)
        self.d2=self.parseD(d2)
        self.l2=self.parseL(l2)
        self.sign=sign
        self.deleted=deleted
        self.coords1=Interval.parseInt(coords1)
        self.coords2=Interval.parseInt(coords2)
    def parseD(self,d):
        if(not rex.find("\[D=(\S+)\]",d)): raise Exception("Can't parse D: "+d)
        return int(rex[1])
    def parseL(self,l):
        if(not rex.find("L=(\S+)",l)): raise Exception("Can't parse L: "+l)
        return int(rex[1])
    def toString(self):
        return self.readID+"\t"+self.guide1.toString()+" [D="+str(self.d1)+\
            "] L="+str(self.l1)+"\t"+self.guide2.toString()+" [D="+\
            str(self.d2)+"] L="+str(self.l2)+"\t"+self.sign+"\t"+self.deleted+\
            "\t"+self.coords1.toString()+"\t"+self.coords2.toString()





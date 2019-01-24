def getLongestSoftMask(rec):
    cigar=rec.getCigar()
    L=cigar.length()
    longest=None
    for i in range(L):
        if(cigar[i].getOp()!="S"): continue
        if(longest is None or cigar[i].getLength()>longest.getLength()):
            longest=cigar[i]
    return longest


def findUnaligned_OLD(unaligned,genome,revGenome):
    L=len(genome)
    querylen=len(unaligned)
    last=L-querylen
    # Search forward strand:
    for i in range(last):
        if(genome[i:(i+querylen)]==unaligned):
            return i
    # Search reverse strand:
    for i in range(last):
        if(revGenome[i:(i+querylen)]==unaligned):
            return L-i-1





    #print(nearestTarget.pos,breakpoint,sep="\t")
    #print(rec.ID,rec.CIGAR.toString(),
          #"mult="+str(rec.flag_hasMultipleSegments()),
          #"aligned="+str(rec.flag_properlyAligned()),
          #"unmapped="+str(rec.flag_unmapped()),
          #"first="+str(rec.flag_firstOfPair()),
          #"second="+str(rec.flag_secondOfPair()),
          #"rev="+str(rec.flag_revComp()),
          #"secondary="+str(rec.flag_secondaryAlignment()),
          #"failed="+str(rec.flag_failedFilters()),
          #"dup="+str(rec.flag_PCRduplicate()),
          #"suppl="+str(rec.flag_supplAlignment()),
          #sep="\t")



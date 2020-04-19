def TrimAlignmentByDiscardingTaxa(alignment,gapTolPerTaxa,treatAmbiguousAsGap=True):
    trimmedAlignment = {}
    numberOfColumns = len(list(alignment.values())[0])
    for taxon in alignment.keys():
        seq = alignment[taxon]        
        if treatAmbiguousAsGap:
            numberOfNonGappedSites = seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C")
        else:
            numberOfGappedSites = seq.count("-")
            numberOfNonGappedSites = numberOfColumns - numberOfGappedSites
        numberOfGappedSites = numberOfColumns - numberOfNonGappedSites
        if (numberOfGappedSites < gapTolPerTaxa*numberOfColumns):
            trimmedAlignment[taxon] = seq
    return (trimmedAlignment)

def TrimAlignmentByDiscardingColumns(alignment,gapTolerancePerColumn,treatAmbiguousAsGap=True):
    trimmedAlignment = {}
    numberOfTaxa = len(list(alignment.keys()))
    numberOfColumns = len(list(alignment.values())[0])
    gapsPerSite = [0] * numberOfColumns
    if treatAmbiguousAsGap:
        for taxon in alignment.keys():
            seq = alignment[taxon]
            for site in range(numberOfColumns):
                if seq[site] == "A" or seq[site] == "T" or seq[site] == "G" or seq[site] == "C":
                    pass
                else:
                    gapsPerSite[site] += 1
    else:
        for taxon in alignment.keys():
            seq = alignment[taxon]
            for site in range(numberOfColumns):
                if seq[site] == "-":
                    gapsPerSite[site] += 1
    for taxon in alignment.keys():
        seq = alignment[taxon]
        trimmedSeq = ""
        for site in range(numberOfColumns):
            if (gapsPerSite[site] < gapTolerancePerColumn * numberOfTaxa):
                trimmedSeq += seq[site]
        trimmedAlignment[taxon] = trimmedSeq        
    return(trimmedAlignment)

def TrimAlignment(alignment,gapTolerancePerColumn,gapTolPerTaxa,treatAmbiguousAsGap=True):
    trimmedAlignmentWithTaxaDiscarded = TrimAlignmentByDiscardingTaxa(alignment,gapTolPerTaxa,treatAmbiguousAsGap=True)
    trimmedAlignment = TrimAlignmentByDiscardingColumns(trimmedAlignmentWithTaxaDiscarded,gapTolerancePerColumn,treatAmbiguousAsGap=True)
    return(trimmedAlignment)
import subprocess as sub

from fileIO import ReadAlignment, WriteAlignment
from alignmentTools import TrimAlignment

partsToRun = {}
partsToRun["parseAlignmentNames"] = True
partsToRun["alignment"] = True
partsToRun["trimAlignment"] = True
partsToRun["phyloInference"] = True
dataPath = "/home/pk/Projects/CoronaTree/data/"
# Concatenate German covid-19 sequences from CHARITE and GISAID, and sequences from Wuhan
seqs_chinese_european_non_german = ReadAlignment(dataPath+"corona_european_china_wholeGenome_trimmed_non_german.fasta")
seqs_german_gisaid = ReadAlignment(dataPath+"gisaid_cov2020_germany_sequences.fasta")
seqs_german_charite = ReadAlignment(dataPath+"charite-SARS-CoV-2.fasta")
seqs_all = seqs_german_gisaid.copy()
seqs_all.update(seqs_german_charite)
seqs_all.update(seqs_chinese_european_non_german)

if partsToRun["parseAlignmentNames"]:
    # replace " " (space), ";" (semicolon), ":" (colon), "," (comma), "()" (parentheses), "'" (quote) with 
    parsedAlignment = {}
    for taxonName in seqs_all.keys():
        seq = seqs_all[taxonName]
        taxonNameWithoutSpecialCharacters = taxonName
        taxonNameWithoutSpecialCharacters = taxonNameWithoutSpecialCharacters.replace(" ", "_")
        taxonNameWithoutSpecialCharacters = taxonNameWithoutSpecialCharacters.replace(";", "_")
        taxonNameWithoutSpecialCharacters = taxonNameWithoutSpecialCharacters.replace(":", "_")
        taxonNameWithoutSpecialCharacters = taxonNameWithoutSpecialCharacters.replace(",", "_")
        taxonNameWithoutSpecialCharacters = taxonNameWithoutSpecialCharacters.replace("(", "_")
        taxonNameWithoutSpecialCharacters = taxonNameWithoutSpecialCharacters.replace(")", "_")
        taxonNameWithoutSpecialCharacters = taxonNameWithoutSpecialCharacters.replace("'", "_")
        parsedAlignment[taxonNameWithoutSpecialCharacters] = seq
    WriteAlignment(parsedAlignment,dataPath+"corona_germany_unaligned.fasta")
if partsToRun["alignment"]:
    command_mafft = "mafft "+dataPath+"corona_germany_unaligned.fasta > " +dataPath+"corona_germany_aligned.fasta"
    sub.call(command_mafft,shell=True)
if partsToRun["trimAlignment"]:    
    alignedSeqs = ReadAlignment(dataPath+"corona_germany_aligned.fasta")
    trimmedAlignment = TrimAlignment(alignedSeqs,gapTolerancePerColumn=0.05,gapTolPerTaxa=0.2,treatAmbiguousAsGap=True)
    WriteAlignment(trimmedAlignment,dataPath+"corona_germany_trimmed.fasta")
if partsToRun["phyloInference"]:
    alignedSeqs = ReadAlignment(dataPath+"corona_germany_trimmed.fasta")
    command_raxml = "/home/pk/Tools/raxml-ng --msa " + dataPath + "corona_germany_trimmed.fasta --tree pars{1} --model GTR --redo"
    sub.call(command_raxml,shell=True)

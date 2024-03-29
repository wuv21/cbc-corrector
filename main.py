from Bio import SeqIO
from Bio import SeqRecord
from typing import Union, Generator
from itertools import combinations,product
from math import prod
import numpy as np
from functools import reduce
from pprint import pprint
from scripts.terminalPrinting import *
import os
import argparse
import pysam
import time
import gzip
import shelve

DNA_ALPHABET = "ATGC"
DNA_EDIT_OPTIONS = {char: [c for c in DNA_ALPHABET if c != char] for char in DNA_ALPHABET}
DNA_EDIT_OPTIONS["N"] = DNA_ALPHABET

def convertPhredToProb(phred: Union[list,str]) -> Union[list,str]:
  """converts phred33 to probability"""
  formula = lambda x: 10 ** (-x / 10)

  if type(phred) is list:
    return list(map(formula, phred))
  
  return formula(phred)
  

def generateSimilarSeqs(seq: str, allowlist: dict, hDist: int = 1) -> list:
  """
  generate nearby neighbors with specified hamning distance
  
  modified from cellranger-atac python code for correcting barcodes
  """ 
  
  canEditIdx = [i for i in range(len(seq)) if seq[i] != 'N']
  mustEditIdx = tuple([i for i in range(len(seq)) if seq[i] == 'N'])

  requiredEdits = len(mustEditIdx)
  similarSeqs = []

  if requiredEdits > hDist:
    return None
  
  # this allows for the condition of x N-containing cbcs with a max hDist of x.
  if requiredEdits == 0:
    minEdits = 1 # setting to one because wouldn't be in this fx if a correction wasn't needed...
  else:
    minEdits = requiredEdits

  for dist in range(minEdits, hDist + 1):
    canEditAllowedN = dist - requiredEdits

    for editIdx in combinations(canEditIdx, canEditAllowedN):
      indices = set(editIdx + mustEditIdx)
      
      edits = product(*["".join(DNA_EDIT_OPTIONS[base]) if i in indices else base for i,base in enumerate(seq)])
      
      for ed in edits:
        editedSeq = "".join(ed)

        if editedSeq in allowlist:
          similarSeqs.append({"edit": editedSeq, "idx": list(indices)})

  return similarSeqs


def main(args):
  tStart = time.time()

  allowlistCbcs = {}
  allowlistCbcsCounter = 0
  uniqAllowlistCbcsCounter = 0
  uniqSusCbcsCounter = 0
  totalCbcsCounter = 0
  susCorrectedCounter = 0

  outputShelve = shelve.open(args.output)
  #finalRecords = []

  with open(args.allowlist) as handle:
    for line in handle:
      line = line.strip()
      allowlistCbcs[line] = 0
  
  susCbcs = {}
  with gzip.open(args.fastq, "rt") as handle:
    for record in SeqIO.parse(handle, "fastq"):
      score = record.letter_annotations["phred_quality"]
      score = convertPhredToProb(score)

      cbc = str(record.seq)
      
      if args.seqWorkflow == "rc":
        cbc = str(record.seq.reverse_complement())
        score.reverse()

      if cbc in allowlistCbcs:
        if allowlistCbcs[cbc] == 0:
          uniqAllowlistCbcsCounter += 1

        allowlistCbcs[cbc] += 1
        allowlistCbcsCounter += 1
        
        # add to finalRecords. The "-1" denotes a valid barcode
        outputShelve[record.id] = cbc + args.BAMTagSuffix
        #finalRecords.append([])

        #finalRecords.append([record.id, cbc + args.BAMTagSuffix])

      else:
        uniqSusCbcsCounter += 1
        #susCbcs[totalCbcsCounter] = {"seq": cbc, "q": score}
        susCbcs[record.id] = {"seq": cbc, "q": score}

        # add to finalRecords. No "-1" denotes a sus barcode.
        #finalRecords.append([record.id, cbc])
        outputShelve[record.id] = cbc

      totalCbcsCounter += 1

  # generate prior distribution of barcodes in allowlist
  allowlistCbcsProb = allowlistCbcs
  for cbc in allowlistCbcs:
    allowlistCbcsProb[cbc] = allowlistCbcs[cbc] / allowlistCbcsCounter

  # find one hamming distance barcodes for each sus barcode
  for scbcRecordIdx in susCbcs:
    scbc = susCbcs[scbcRecordIdx]["seq"]
    q = susCbcs[scbcRecordIdx]["q"]

    potentialEdits = generateSimilarSeqs(seq = scbc, allowlist = allowlistCbcs, hDist = 1)
    
    if potentialEdits == None:
      continue

    # generate posterior probs
    posProbs = []
    posEdits = []
    for editedScbc in potentialEdits:
      
      correctBc = editedScbc["edit"]
      seqError = prod([q[i] for i in editedScbc['idx']])
      pp = allowlistCbcsProb[correctBc] * seqError
      
      posProbs.append(pp)
      posEdits.append(editedScbc)
      

    if len(posProbs) > 0 and sum(posProbs) > 0:
      posProbs = np.array(posProbs)
      posProbs = posProbs / np.sum(posProbs)

      # get the max edit that is > 0.975. if there's a tie, barcode does not get corrected.
      maxPosProb = np.amax(posProbs)
      if maxPosProb < args.posProbThresh:
        continue

      # debugging
      # print('--------')
      # print(scbc)
      # print(potentialEdits)
      # print([allowlistCbcsProb[cbc['edit']] for cbc in potentialEdits])
      # print(correctBc)
      # print(posProbs)

      maxPosProbIdxs = [i for i,x in enumerate(posProbs) if x == maxPosProb]
      if len(maxPosProbIdxs) == 1:
        susCorrectedCounter += 1

        outputShelve[scbcRecordIdx] = posEdits[maxPosProbIdxs[0]]["edit"] + args.BAMTagSuffix
        # finalRecords[scbcRecordIdx][1] = posEdits[maxPosProbIdxs[0]]["edit"] + args.BAMTagSuffix

  # output list
  #with open(args.output, "w") as outfile:
  #  for i in range(0, len(finalRecords)):
  #    outfile.write("\t".join(finalRecords[i]))
  #    outfile.write("\n")

  print("Finished correcting. Adding tag to BAM file")
  
  # convert finalRecords to dictionary
  #recordDict = {rec[0]: rec[1] for rec in finalRecords}

  # process bam file
  origBam = pysam.AlignmentFile(args.bam, "rb")
  outBam = pysam.AlignmentFile(args.outBam, "wb", template = origBam)

  for read in origBam.fetch(until_eof = True):
    qname = read.query_name
    if qname in outputShelve:
      read.set_tag(args.BAMTagField, outputShelve[qname], value_type = "Z")
    
    outBam.write(read)

  origBam.close()
  outBam.close()
  outputShelve.sync()
  outputShelve.close()

  print("Time elapsed: {} seconds".format(str(time.time() - tStart)))
  print("Unique barcodes found in allowlist: {}".format(uniqAllowlistCbcsCounter))
  print("Total barcodes found in allowlist: {}".format(allowlistCbcsCounter))
  print("Total barcodes that were sus: {}".format(uniqSusCbcsCounter))
  print("Total sus barcodes that were corrected: {}".format(susCorrectedCounter))
  print("Total barcodes processed: {}".format(totalCbcsCounter))
  print()


if __name__ == "__main__":
  # set up command line arguments
  parser = argparse.ArgumentParser(
    description = "Extract, compile, and correct cell barcodes from single-cell applications using 10X barcode correction algorithm")

  parser.add_argument("--fastq",
    required = True,
    help = "Fastq.gz file containing barcode reads")
  parser.add_argument("--allowlist",
    required = True,
    help = "Txt file containing allowlist barcodes (one per line)")
  parser.add_argument("--output",
    required = True,
    help = "Output file path")
  parser.add_argument("--bam",
    required = True,
    help = "BAM file to assign corrected barcodes")
  parser.add_argument("--outBam",
    required = True,
    help = "BAM output file with corrected barcodes")
  parser.add_argument("--BAMTagField",
    default = "CB",
    help = "Prefix to add to BAM cell barcode tag")
  parser.add_argument("--BAMTagSuffix",
    default = "-1",
    help = "Suffix to add to BAM cell barcode tag")
  parser.add_argument("--posProbThresh",
    default = 0.975,
    type = float,
    help = "Minimum probability needed for barcode correction")
  parser.add_argument("--seqWorkflow",
    default = "rc",
    choices = ["rc", "fwd"],
    help = "Illumina sequencing workflow. Default is the reverse complement workflow (i.e. NovaSeq v1.5 kits)")

  args = parser.parse_args()
  
  if not os.path.exists(args.fastq):
    raise Exception("Fastq file not found")

  if not os.path.exists(args.allowlist):
    raise Exception("Allowlist text file not found")
  
  if not os.path.exists(os.path.dirname(args.output)):
    raise Exception("Directory for output text file not found")

  if not os.path.exists(os.path.dirname(args.outBam)):
    raise Exception("Directory for output bam file not found")

  main(args)

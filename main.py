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

DNA_ALPHABET = "ATGC"
DNA_EDIT_OPTIONS = {char: [c for c in DNA_ALPHABET if c != char] for char in DNA_ALPHABET}
DNA_EDIT_OPTIONS["N"] = DNA_ALPHABET

def convertPhredToProb(phred: Union[list,str]) -> Union[list,str]:
  """converts phred33 to probability"""
  formula = lambda x: 10 ** (-x / 10)

  if type(phred) is list:
    return list(map(formula, phred))
  
  return formula(phred)
  

def generateSimilarSeqs(seq: str, allowlist: dict, hDist: int = 1) -> Generator[dict, None, None]:
  """
  generate nearby neighbors with specified hamning distance
  
  modified from cellranger-atac python code for correcting barcodes
  """ 
  
  canEditIdx = [i for i in range(len(seq)) if seq[i] != 'N']
  mustEditIdx = tuple([i for i in range(len(seq)) if seq[i] == 'N'])

  requiredEdits = len(mustEditIdx)
  if requiredEdits > hDist:
    return None

  for dist in range(requiredEdits + 1, hDist + 1):
    for editIdx in combinations(canEditIdx, dist - requiredEdits):
      indices = set(editIdx + mustEditIdx)
      
      
      edits = product(*["".join(DNA_EDIT_OPTIONS[base]) if i in indices else base for i,base in enumerate(seq)])
      
      for ed in edits:
        editedSeq = "".join(ed)

        if editedSeq in allowlist:
          yield {"edit": editedSeq, "idx": list(indices)}


def main(args):
  allowlistCbcs = {}
  allowlistCbcsCounter = 0
  uniqAllowlistCbcsCounter = 0
  uniqSusCbcsCounter = 0
  totalCbcsCounter = 0
  susCorrectedCounter = 0


  finalRecords = []

  with open(args.allowlist) as handle:
    for line in handle:
      line = line.strip()
      allowlistCbcs[line] = 0
  
  susCbcs = {}
  with open(args.fastq) as handle:
    for record in SeqIO.parse(handle, "fastq"):
      score = record.letter_annotations["phred_quality"]
      score = convertPhredToProb(score)

      cbc = str(record.seq)

      # double check if cbc is rc, if so also need to reverse order of q scores
      #cbc = str(record.seq.reverse_complement())
      #score.reverse()

      if cbc in allowlistCbcs:
        if allowlistCbcs[cbc] == 0:
          uniqAllowlistCbcsCounter += 1

        allowlistCbcs[cbc] += 1
        allowlistCbcsCounter += 1
        
        # add to finalRecords. The "-1" denotes a valid barcode
        finalRecords.append([record.id, cbc + "-1"])

      else:
        uniqSusCbcsCounter += 1
        susCbcs[totalCbcsCounter] = {"seq": cbc, "q": score}

        # add to finalRecords. No "-1" denotes a sus barcode.
        finalRecords.append([record.id, cbc])

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

      maxPosProbIdxs = [i for i,x in enumerate(posProbs) if x == maxPosProb]
      if len(maxPosProbIdxs) == 1:
        susCorrectedCounter += 1
        finalRecords[scbcRecordIdx][1] = posEdits[maxPosProbIdxs[0]]["edit"]

  # output list
  with open(args.output, "w") as outfile:
    for i in range(0, len(finalRecords)):
      finalRecords[i] = "\t".join(finalRecords[i])
    
    outfile.write("\n".join(finalRecords))
  
  print("Unique barcodes found in allowlist: {}".format(uniqAllowlistCbcsCounter))
  print("Total barcodes found in allowlist: {}".format(allowlistCbcsCounter))
  print("Total barcodes that were sus: {}".format(uniqSusCbcsCounter))
  print("Total sus barcodes that were corrected: {}".format(susCorrectedCounter))
  print("Total barcodes processed: {}".format(totalCbcsCounter))
  

if __name__ == "__main__":
  # set up command line arguments
  parser = argparse.ArgumentParser(
    description = "Extract, compile, and correct cell barcodes from single-cell applications using 10X barcode correction algorithm")

  parser.add_argument("--fastq",
    required = True,
    help = "Fastq file containing barcode reads")
  parser.add_argument("--allowlist",
    required = True,
    help = "Txt file containing allowlist barcodes (one per line)")
  parser.add_argument("--output",
    required = True,
    help = "Output file path")
  parser.add_argument("--posProbThresh",
    default = 0.975,
    type = float,
    help = "Minimum probability needed for barcode correction")

  args = parser.parse_args()
  
  if not os.path.exists(args.fastq):
    raise Exception("Fastq file not found")

  if not os.path.exists(allowlist):
    raise Exception("Allowlist text file not found")

  main(args)

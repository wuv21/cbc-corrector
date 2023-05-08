# cbc-corrector

## About
cbc-corrector is a Python tool to extract cell barcodes from a fastq file in order to add to alignment records. For 10X Genomics applications, an allowlist of cell barcodes (i.e. these barcodes are expected) allows for end users to identity cell barcodes that may be otherwise excluded due to a sequencing error. This tool identifies cell barcodes that are 1 Hamming distance away from a cell barcode in the allowlist. This tool uses the Bayesian strategy and algorithm described by 10X Genomics (see [here](https://kb.10xgenomics.com/hc/en-us/articles/115003822406-How-does-Cell-Ranger-correct-barcode-sequencing-errors-)) and is derived from the code found in the [cellranger-atac pipeline](https://github.com/10XGenomics/cellranger-atac). 

## Getting started

### Requirements
1. See `env.yaml` for the conda environment and prerequesites.
2. Create a conda env using the yaml file: `conda env create -f env.yml -n cbc-correction`
3. Download allowlist from 10X Genomics or other sources.
4. Run cbc-corrector (see below)

```bash
# run cbc-corrector
python main.py \
  --fastq=tmp_R2.fastq.gz \
  --allowlist=737K-cratac-v1.txt \
  --output=tmp_barcodes.txt \
  --bam=tmp.bam \
  --outBam=tmp.corrected.bam
```

## Parameters

- `--fastq` *(required)* Fastq.gz file containing barcode reads. For 10X Genomics' ATAC assay, this is typically the R2 fastq file (note the different nomenclature due to Illumina's bcl2fastq naming system: R1 = read 1, I1 = index 1, R2 = index 2, R3 = read 2). Note that this file needs to be gunzipped/compressed.
- `--allowlist` *(required)* Txt file containing allowlist barcodes (one per line). For scATAC assays, the file can be found in the cellranger-atac package and is named "737K-cratac-v1.txt".
- `--output` *(required)* Output file path.
- `--posProbThresh` Default = 0.975. Float value that denotes the minimum probability needed for barcode correction.
- `--bam` *(required)* Input bam file that needs assignment of corrected cell barcode.
- `--outBam` *(required)* Output bam file path with corrected cell barcodes.
- `--BAMTagField` Default = "CB". Tag field to write corrected cell barcode to in output BAM file.
- `--BAMTagSuffix` Default = "-1". Suffix to write to corrected cell barcode to denote correct cell barcode. This default matches with the 10X Genomics nomenclature.
- `--seqWorkflow` Default = "rc". Choose either "rc" or "fwd" to match the sequencing workflow used. NovaSeq v1.5 kits used the "rc" (reverse complement) workflow while v1 kits use the "fwd" (forward strand) workflow. 


## Outputs
- Text file that contains two tab-separated columns. First column is the read name. Second column is the extracted barcode. If the barcode has a "-1" suffix, this means that the barcode is in the allowlist OR is the corrected barcode (also is found in the allowlist).

See below for an example
```
A00626:220:HL2KFDMXX:1:1101:19180:1047	NAGGAGCAGACCATAA
A00626:220:HL2KFDMXX:1:1101:20464:1047	NCACTCGTCGCTCTAC
A00626:220:HL2KFDMXX:1:1101:21856:1047	NCAGAAAAGCGAGCTA
A00626:220:HL2KFDMXX:1:1101:23032:1047	NTGTTCGAGTTAGCAA
A00626:220:HL2KFDMXX:1:1101:24044:1047	NCCCAGAGAGATTACA
A00626:220:HL2KFDMXX:1:1101:24912:1047	NAGCATGAGGCTGGAT
A00626:220:HL2KFDMXX:1:1101:26268:1047	NAATGAGGTGGACTGA
A00626:220:HL2KFDMXX:1:1101:26485:1047	NGCGTAACAACGAGGT
A00626:220:HL2KFDMXX:1:1101:26648:1047	NTATGTGATTAGAACG
A00626:220:HL2KFDMXX:1:1101:27751:1047	NGAGTCATCAGGATCT
A00626:220:HL2KFDMXX:1:1101:27896:1047	NCTGTCCTCGGACGAA
A00626:220:HL2KFDMXX:1:1101:28167:1047	NGCCGCAAGCCATTCA
A00626:220:HL2KFDMXX:1:1101:28800:1047	NAGGTCCTCCATAACG
A00626:220:HL2KFDMXX:1:1101:29631:1047	NTCTAACTCAGATACC
A00626:220:HL2KFDMXX:1:1101:30463:1047	NAATCTGCAGACTAAA
A00626:220:HL2KFDMXX:1:1101:32090:1047	NCTTGCAAGAATCAAC
A00626:220:HL2KFDMXX:1:1101:32163:1047	NACATCCAAGGCGTAT
A00626:220:HL2KFDMXX:1:1101:32325:1047	NCATTGATCCATGACA
A00626:220:HL2KFDMXX:1:1101:1018:1063	GCGCCAAGTACAAGCG-1
A00626:220:HL2KFDMXX:1:1101:1488:1063	CGAGTTAGTGCATCAT-1
A00626:220:HL2KFDMXX:1:1101:1976:1063	TTGCGGGAGAAATACC-1
A00626:220:HL2KFDMXX:1:1101:2013:1063	ACCAAACTCAGTACAC-1
A00626:220:HL2KFDMXX:1:1101:2790:1063	TTACCGCCATGGGACA-1
A00626:220:HL2KFDMXX:1:1101:3043:1063	AGATTCGTCCAATAGC-1
A00626:220:HL2KFDMXX:1:1101:3423:1063	GCAACCGTCTACCCGT-1
A00626:220:HL2KFDMXX:1:1101:4273:1063	CGTGGCACAAGGATGC-1
A00626:220:HL2KFDMXX:1:1101:5249:1063	TACGCCTGTCAACGGA-1
A00626:220:HL2KFDMXX:1:1101:5303:1063	CAGTGCGTCGCAGATT-1
A00626:220:HL2KFDMXX:1:1101:5466:1063	GAGTGAGTCAAGAGGC-1
```

- BAM file with corrected cell barcodes inserted for each record if it contains a valid cell barcode.

# synthetic_reads
Generate synthetic Illumina reads using a provided multifasta.

## Requirements

tested with 
- python 3.11.0
- numpy 1.23.5

## Usage

Basic usage requires only the input of a fasta sequence file

`python3 synthetic reads.py -f path/to/sequence.fasta`

## Procedure

1. choose a random sequence in the fasta (called a contig herein)
2. choose a random position in the chosen contig
3. select a DNA fragment size by choosing a random number between the minimum and maximum fragment length (set using `-i`)
4. construct two read sequences as follows (read length set using `-l`):
   - if both ends of the dna fragment are within the selected contig, read sequences are just the forward and reverse complement sequence of the contig between the two points
   - if one or both reads extend off the end of the contig by fewer than the read length - 50 bases (i.e., more than 50 bases of the read will be in the contig), use contig sequence and make it up to the read length bases with random sequence
   - if one or both reads are off the end of the contig by more than the read length - 50 bases, generate a random sequence
5. with 50% probability, reverse which read is going to be written to which file so we don’t have one file of all forward and one of all reverse reads when mapped
6. generate random quality scores by sampling from a normal distribution with user controlled mean and sd (set with `-q` and `-s`, respectively)
7. write the reads to two files. Filenames are made by taking the user-provided output prefix (`-o`) and appending either “1.fastq” or “2.fastq”

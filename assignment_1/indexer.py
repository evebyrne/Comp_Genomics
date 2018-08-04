#Assignment 1 indexer.py
#Student Num 929395
#April 7th 2018
##!/usr/bin/env python
import pysam
import sys

outfile = open("index.txt", "w")

#Check the command line arguments
parameters = sys.argv
argc = len(parameters)
if argc != 3:
  print "Usage: <reference file (fasta)> <k>"
  sys.exit(0)

reference = ''
f = pysam.FastxFile(sys.argv[1])
for ref in f:
    reference = ref.sequence.upper()

k = int(sys.argv[2])
index = {}
#print [reference[i:i+k] for i in range(len(reference)-k+1)]
for i in range(len(reference)-k+1):
    kmer = reference[i:i+k]
    if kmer not in index:
         index[kmer] = [i]
    else: index[kmer].append(i)

outfile.write("INDEX: "+ref.name+"\n")
keylist = index.keys()
keylist.sort()
for key in keylist:
    outfile.write("%s %s\n" % (key, " ".join(str(x) for x in index[key])))
outfile.close()

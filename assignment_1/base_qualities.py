#!/usr/bin/env python
#Assignment 1
#Student Num 929395
#April 7th 2018

from Bio import SeqIO
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pysam

#Support class to execute the handiwork of read alignment to a reference sequence
class Aligner1:
  def __init__(self, ref):
    self.refname = ref.id
    self.seq = ref.seq.upper()
    self.rseq = ref.seq.reverse_complement()
    self.length = len(ref.seq)
  def align(self, read):
    alignments = []
    return alignments

# compute hamming distance equal length strings
def hamming_distance(a, b):
    assert(len(a) == len(b))
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))

def qualities(read, ref, strand, quality):
    qs = read.letter_annotations["phred_quality"]
    assert(len(read.seq)==len(ref))
    if (strand == '-'):
        zz = zip(read.seq.reverse_complement(), ref, qs[::-1])
    else : zz =zip(read.seq, ref, qs)
    #print "LEN " + str(len(zz))
    j = 0
    for rd, rf, q in zz:
        if rd != rf:
           if strand == '-':
               if (len(read.seq)-j-1+1) not in quality:
                     #print "- j: %s, len-j-1: %s " %(j , len(read.seq)-j-1)
                     quality[len(read.seq)-j-1+1] = [q]
               else: quality[len(read.seq)-j-1+1].append(q)
           else : 
               if j+1 not in quality:
                     #print "+ j: " + str(j)
                     quality[j+1] = [q]
               else: quality[j+1].append(q)
           #print str(rd) + " " + str(rf) + " " + str(q)
           #print quality
        j +=1       

#Check the command line arguments
parameters = sys.argv
argc = len(parameters)
if argc != 4:
  print "Usage: <reference file (fasta)> <read file (fasta)> <max hamming distance>"
  sys.exit(0)

#Read the reference sequence and initiate the aligner
try:
  for s in SeqIO.parse(sys.argv[1], "fasta"):
    ref = Aligner1(s)
    #print ref.seq
    refname = ref.refname
    break #Stop after the first sequence in the reference
except IOError as e:
  print "Could not read reference sequence file (see below)!"
  print e
  sys.exit(1)


#open and parse the index file
inputfile = open('index.txt', 'r')
index = {}
first_line = inputfile.readline()
for line in inputfile:
     line = line.strip()
     l = line.split()
     kmer = l[0]
     rest = l[1:]
     index[kmer] = map(int, rest)
     #print index[kmer]
k = len(kmer)

outfile = open('alignment.txt', 'w')
outfile.write('READ_NAME \t REF_NAME \t POS \t STRAND \t NUMBER_OF_ALIGNMENTS \t HAMMING_DISTANCE \n')

try:
        qs = []
        qualitymismatches = {}
        alignment_ham = "*"
#READ LOOP
	for r in SeqIO.parse(sys.argv[2],"fastq"):
                qs.append(r.letter_annotations["phred_quality"])
		hamming = sys.argv[3]
		pos = 0
		count = 0
		best_pos = -1
		strand = '+'
		visited = []
		rvisited = []
		read = Aligner1(r)

		#KMER LOOP
		for i in range(read.length-k+1):
		       kmer = read.seq[i:i+k] 
		       rkmer = read.rseq[i:i+k]
		       if kmer in index:
		          #POSITIONS LOOP - NORMAL
		          for pos in index[kmer]:
		              inref = ref.seq[int(pos)-i:(int(pos)+k+(read.length-k)-i)]
		              if (int(pos)-i) not in visited:
		                 if read.length == len(inref):
		                      diff1 = hamming_distance(read.seq, inref)
		                      if int(diff1) <= int(hamming):
		                           if(int(diff1)==int(hamming) and best_pos != -1):
		                                  count +=1
		                                  if ((pos-i) < best_pos):
		                                        best_pos = int(pos)-i	  
		                           else:  
		                                  best_pos = int(pos)-i
		                                  count =1
		                                  hamming = diff1   
		                 visited.append(int(pos)-i)
		       if rkmer in index:
		          #POSITIONS LOOP - REVERSECOMP
		          for pos in index[rkmer]:
		              inref = ref.seq[int(pos)-i:(int(pos)+k+(read.length-k)-i)]
		              if (int(pos)-i) not in rvisited:
		                 if read.length == len(inref):
		                      diff2 = hamming_distance(read.rseq, inref)
		                      if int(diff2) <= int(hamming):
		                           if(int(diff2)==int(hamming) and best_pos != -1):
		                                  count +=1
		                                  if ((pos-i) < best_pos):
		                                        best_pos = int(pos)-i
		                                        strand = '-'
		                           else:  
		                                  best_pos = int(pos)-i
		                                  count =1
		                                  hamming = diff2
		                                  strand = '-'
		                 rvisited.append(int(pos)-i)
                #PALINDROMES
                if (read.seq==read.rseq and count >0):
                       count+=1
                alignment_ham = str(hamming)
                if count == 0:
                     alignment_ham = strand ='*'
                
                if (hamming > 0 and best_pos != -1):
                     qualities(r, ref.seq[int(best_pos):int(best_pos)+read.length], strand, qualitymismatches)
                alignment_pos = best_pos + 1
		
                
		output = read.refname + "\t" + refname + "\t" + str(alignment_pos) + "\t" + str(strand) + "\t" + str(count) + "\t" + alignment_ham + "\n"
		outfile.write(output)
        #print str(qualitymismatches)
        qmm = []
        for key in qualitymismatches.keys():
            #print str(qualitymismatches[key])
            qmm.append(qualitymismatches[key])
        plt.figure(0)
        #print [list(i) for i in zip(*qs)]
        #print qmm
        plt.boxplot([list(i) for i in zip(*qs)])
        plt.savefig('QS.png')
        plt.figure(1)
        plt.boxplot(qmm)
        plt.savefig('QSMM.png')
except IOError as e:
    print 'error!'
    print e
    sys.exit(1)
inputfile.close()


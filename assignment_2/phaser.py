import pysam
import sys

#BAM = 0 based
#reads from vcf file are coming from the pilupcolumn.pos method which comes from bam so vcf pos are 0 based
bamfile = pysam.AlignmentFile(sys.argv[1],'rb')
f = pysam.VariantFile('task1.vcf')
prevpos = 0
het = []

for line in f:
     if not line:
        flag = false
     info = line.info.keys()[0]
     if float(info) < float(0.8):
         het.append(line)
         #print het
i = 0
bamfile = pysam.AlignmentFile(sys.argv[1],'rb')  
haplotypes = {}
phased = False      
for line in het:
       phased = False
       r1 = bamfile.fetch("chr15", line.pos, line.pos+1)
       rn1 = set([read for read in r1])
       j = i+1  
       while not phased and j < len(het):
         line2 = het[j]
         if line2.pos > line.pos: 
	       r2 = bamfile.fetch("chr15", line2.pos, line2.pos+1)
	       rn2 = set([read for read in r2])
	       if len(rn1.intersection(rn2)):
		    haplotypes = {}
		    inter = rn1.intersection(rn2)
		    for read in inter:
			 base1 = ""
			 base2 = ""
			 pos1 = line.pos - read.pos
			 pos2 = line2.pos - read.pos
			 if len(read.cigar) != 1:
			      maxpos = 0
			      cigar = read.cigar
			      lencigar = len(cigar)
			      ind = 0
			      if (cigar[0])[0] == 4:
			           lengthsoft = (cigar[0])[1]
			           pos1 += lengthsoft
			           pos2 += lengthsoft
			           ind +=1
			      while int(ind) < int(lencigar):
			           code = (cigar[ind])[0]
			           value = (cigar[ind])[1]
			           if code == 0:
			                maxpos = maxpos + value
			           elif code == 1:
			                  pos1 += value
			                  pos2 += value
			           elif code == 2:
			                  pos1 -= value
			                  pos2 -= value 
			           if pos1 <= maxpos and str(base1) == str(""):
			                base1 = read.seq[pos1]
			           if pos2 <= maxpos and str(base2) == str(""): 
			                base2 = read.seq[pos2]
			           ind += 1
			 #normal CIGAR                                                                       
			 else: 
			     base1 = read.seq[pos1]
		      
			     if int(pos2) < int(len(read.seq)):
			        base2 = read.seq[pos2]
			 haplotype = base1 + base2
			 if haplotype not in haplotypes: 
			        haplotypes[haplotype] =1
			 else: haplotypes[haplotype] +=1
		    h1 = line.alleles[0]+line2.alleles[0]
		    h2 = line.alleles[0]+line2.alleles[1]
		    h3 = line.alleles[1]+line2.alleles[0]
		    h4 = line.alleles[1]+line2.alleles[1]
		    total = 0
		    prob = {}
		    p1 = p2 = 0
		    phased = False
		    for key in haplotypes:
		       if key == h1 or key == h2 or key == h3 or key == h4:
			  total += haplotypes[key]
		    for key in haplotypes:
		       if key == h1 or key == h2 or key == h3 or key == h4:
			  prob[key] = float(haplotypes[key])/total
		    for key in prob:
		       if prob[key] >= 0.9:
			  phased = True
		       if key == h1 or key == h4:
			  p1 += prob[key]
		       else: p2 += prob[key]
		    if format(float(p1), '.2f') >=format(float(0.9), '.2f'):
	                  haplotype = str(h1) + "/" + str(h4)
			  phased = True
		    elif format(float(p2), '.2f') >=format(float(0.9), '.2f'):
			  haplotype = str(h2) + "/" + str(h3)
			  phased = True 
		    else: 
			  haplotype = "X"
                          phased = False
         j +=1
       if phased:
          print "pos: " + str(line.pos) + " " + str(line2.pos)+" alleles: " + str(line.alleles) +" "+ str(line2.alleles)+"\npossible haplotypes:"+str(prob.keys())+"\n"+"prob: " + str(prob) +"\n"+"haplotype: " + str(haplotype)+"\n"
       i +=1





















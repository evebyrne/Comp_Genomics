import pysam
import sys

#positions in fastq are 1 based
f = pysam.FastxFile('reference_human.fa')
for read in f: 
    ref = (read.sequence).upper()
    ref_pos = read.name
pos =  (ref_pos.split(":")[1]).split('-')
start = pos[0]
end = pos[1]
outfile = open('task1.vcf', 'w')
outfile.write("##fileformat=VCFv4.0\n")
outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
bamfile = pysam.AlignmentFile(sys.argv[1],'rb')
keys = ['A', 'C', 'G', 'T', 'N']

for pileupcolumn in bamfile.pileup():
    read_count = pileupcolumn.n
    #positions in BAM are 0 based
    ref_pos = pileupcolumn.pos +1
    ref_index = ref_pos - int(start)
    ref_base = -1
    if int(ref_pos) >= int(start) and int(ref_pos) <= int(end): 
        ref_base = ref[ref_index]
    base_dict = {'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
    #quality score
    score_dict = {'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
    prob_dict = {'A':[0,0], 'T':[0,0], 'C':[0,0], 'G':[0,0], 'N':[0,0], 'X':[]}
    
    count = read_count

    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            #count the observed nucleotides
            #count the quality scores
            score = pileupread.alignment.query_qualities[pileupread.query_position]
            if score >= 20: 
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                base_dict[base] += 1
                score_dict[base] += score
                
            else: 
              count = count -1 
              base = pileupread.alignment.query_sequence[pileupread.query_position]
#sometimes the probablity < 1 even if all same bases - > think this  has to do with the if not is_del or refskip thing, as in some reads are not being counted but are included in the read_count which will mean / by a larger num than actually counted so prob < 1 even if all bases the same 
# this is probably the reason for the if base_count >=2 which has to make sure that the probability is actually between 0.2 and 0.8 because there is >1 base rather than some reads(bases really) #not being counted 
    base_count = 0
    genotype = ""
    flag = 0
    for key in keys:
        
        if base_dict.get(key) > 0:
            prob = float(base_dict.get(key))/count
            
            if prob >= 0.2 and str(key) != str(ref_base) and ref_base != -1:
              prob_dict[key][0] = prob
              prob_dict[key][1] = score_dict[key]/float(base_dict[key])
              genotype = str(key) + str(key)
              quality = prob_dict[key][1] 
              alt = key
              AF = prob
              var = key
              if prob <=0.80:
                 flag = 1 
    if flag == 1:
       for key in keys: 
            if base_dict.get(key) > 0:
               prob = float(base_dict.get(key))/count
               if prob >=0.20 and prob <= 0.80 and str(key) != str(var):
                      genotype = str(var) + str(key)         
            
    if genotype != "":
       output = str(pileupcolumn.reference_name)+'\t'+str(pileupcolumn.pos)+'\t.\t'+str(ref_base)+'\t'+str(alt)+'\t'+str(prob_dict[var][1])+'\tPASS\t'+str(prob_dict[var][0])+'\tPASS\t'+str(genotype)+'\n'
       outfile.write(output)

outfile.close()

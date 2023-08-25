#!/usr/bin/env python3
import string
import random



def random_exon_id(length=20):
    input_chars = list(string.digits) + list(string.ascii_uppercase)
    return 'RAND' + ''.join(random.choice(input_chars) for i in range(length))

random_exon_id(10)

FEATURE=2
ANNOTATION=8
hsa_filename = 'star_index/hsa_hairline.gtf'
with open('star_index/hsa_hairline_gene_update.gtf', 'w') as ofile:
    for line in open(hsa_filename):
        gtf_line = line.split('\t')
        if gtf_line[FEATURE] =='transcript':
            gene_line = line.split(";")[0] + "\n"
            gene_line = gtf_line.copy()
            gene_line[ANNOTATION] = gene_line[ANNOTATION].split(';')[1].strip()
            gene_line[FEATURE] = 'gene'
            gene_line = '\t'.join(gene_line) + '\n'
            
            transcript_line = line
            exon_line = gtf_line.copy()
            exon_line[ANNOTATION] = f'exon_id "{random_exon_id()}"; ' + exon_line[ANNOTATION]
            exon_line[FEATURE] = 'exon'
            exon_line = '\t'.join(exon_line)
            ofile.write(gene_line)
            ofile.write(transcript_line)
            ofile.write(exon_line)
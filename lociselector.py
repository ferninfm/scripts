#! /usr/bin/env python
#----------------------------------------------#
# Script to survey radtag catalog              #
#                                              #
# Author: Emiliano Trucchi, modified FFM       #
#                                              #
#----------------------------------------------#
# A. Define loci_selector function: FFM added snp_pos_min and GC cont min
#
def loci_selector(input_file,snp_pos_min,snp_pos_max,het_max_allowed,GC_cont_min,GC_cont_max,max_missing_ind):
#
	file = open(input_file, "r")
#1 get the header
	header_raw = file.readline()
	header = header_raw.split("\n")
	header = header[0].split("\t")
#2 initialize output file
	output_selection = input_file
	output_selection += "_selected"
	out = open(output_selection, "w")
	out.write(header_raw)
#3 set the counters FFM
	counter_global_loci = 0
	counter_too_many_missing_drop = 0
	counter_3all_drop = 0
	counter_loci_dropped = 0
	counter_het_dropped = 0
	counter_deleveraged = 0
	counter_invariant_drop = 0
	counter_GChigh = 0
    counter_GClow = 0 #FER
	counter_Ns_dropped = 0
	invariant = 0
#4 set the returned object
	global_position_list = []
	cat_id_list = []
	dropped_loci = []
	list_sample_size = []
	GC_content_list = []
	list_ind_by_locus = []

#5 start the loop on the lines with data
	for line in file.readlines():
		if line == '\n':
			break
		else:
			counter_global_loci += 1
			row = line.split("\n")
			row = row[0].split("\t")
#6 get field values
			cat_id = row[0]
			annotation = row[1]
			chrom = row[2]
			bp = row[3]
			consensus = row[4]
			num_par = row[5]
			num_prog= row[6]	
			num_snp = int(row[7]) # Used in 8
			snp_position = row[8] ######## Watch out
			num_all = row[9]
			alt_alleles = row[10]
			deleveraged = row[11]
			geno = row[12:len(row)]
			hetero = ''.join(geno)
			Cs = Decimal(consensus.count('C'))
			Gs = Decimal(consensus.count('G'))
			GC = (Cs+Gs)/len(consensus)
			
#7 get GC_content of the locus
			Cs = Decimal(consensus.count('C'))
			Gs = Decimal(consensus.count('G'))
			GC = (Cs+Gs)/len(consensus)

#8 check and remove SNPs above a fixed position and correct snp_position, num_snp and alleles field
			if num_snp < 1:
				invariant += 1
			else:
				snp_position = snp_position.split(";")
				position_index = []
				i = 0
				for snp in snp_position:
					snp = snp.split(",")
					position = int(snp[0])
					## CHANGED FFM if position < snp_pos_max:
                    if snp_pos_min < position < snp_pos_max:
						position_index.append(i)
						global_position_list.append(position) 
					i += 1
				snp_position_corr = ";".join(snp_position[0:len(position_index)])
				num_snp_corr = len(position_index)

				
				alleles = alt_alleles.split(";")
				new_allele_list = []
				for allele in alleles:
					new_allele = allele[0:len(position_index)]
					new_allele_list.append(new_allele)
				alleles_corr = ";".join(set(new_allele_list))
				num_all_corr = len(set(new_allele_list))

#9 remove loci with number of missing data above threshold
				missing_ind_by_locus = geno.count("")
				real_ind_num_by_locus = len(geno)-geno.count("")
				list_ind_by_locus.append(real_ind_num_by_locus)
				if real_ind_num_by_locus < max_missing_ind:
					counter_loci_dropped += 1
					dropped = row[0]
					dropped_loci.append(dropped)
					counter_too_many_missing_drop += 1					

#10 remove loci with no SNP before the max snp position allowed
				elif num_snp_corr == 0:
					counter_loci_dropped += 1
					dropped = row[0]
					dropped_loci.append(dropped)
					counter_invariant_drop += 1	

#11 remove loci with GC content over a given threshold
				elif GC > float(GC_cont_max):
					counter_loci_dropped += 1
					dropped = row[0]
					dropped_loci.append(dropped)
					counter_GChigh += 1
					

#12 check and remove loci with more than two alleles in any of the samples
				elif any(len(ind) > num_snp*2+1 for ind in geno):
					counter_loci_dropped += 1
					dropped = row[0]
					dropped_loci.append(dropped)
					counter_3all_drop += 1

#13 check and remove loci deleveraged in stacks
				elif int(deleveraged) != 0:
					counter_loci_dropped += 1
					dropped = row[0]
					dropped_loci.append(dropped)
					counter_deleveraged += 1

#14 check and remove loci with heterozygosity higher than selected
				elif hetero.count('/') > (len(geno)-geno.count(''))*het_max_allowed: 
					counter_loci_dropped += 1
					dropped = row[0]
					dropped_loci.append(dropped)
					counter_het_dropped += 1

#15 remove loci with Ns in any of the individual
				elif ''.join(alleles_corr.split(';')).count('N') != 0:
					counter_loci_dropped += 1
					dropped = row[0]
					dropped_loci.append(dropped)
					counter_Ns_dropped += 1
#16 write the line of selected locus in the output file and close input and output
				else:
					cat_id_list.append(cat_id)
					GC_content_list.append(GC)
					
					sample_size = len(filter(bool, geno))
					list_sample_size.append(sample_size)

					out.write(cat_id+'\t'+annotation+'\t'+chrom+'\t'+bp+'\t'+consensus+'\t'+num_par+'\t'+num_prog+'\t'+str(num_snp_corr)+'\t'+snp_position_corr+'\t'+str(num_all_corr)+'\t'+alleles_corr+'\t'+deleveraged+'\t')
					
					j = 1
					for ind in geno:
						ind = ind.split('/')
						new_elem_list = []
						for elem in ind:
							new_elem = elem[0:len(position_index)]
							new_elem_list.append(new_elem)
						real_alleles = set(new_elem_list)
						ind_corr = "/".join(real_alleles)
						if j < len(geno):
							out.write(ind_corr+'\t')
						elif j == len(geno):
							out.write(ind_corr+'\n')
						j += 1 
					

	file.close()
	out.close()
	print "Number of analyzed loci:", counter_global_loci
	print "Number of invariant loci:", invariant
	print "Number of loci with less than", max_missing_ind, "individuals:", counter_too_many_missing_drop
	print "Number of loci without polymorphisms after removing terminal SNPs:", counter_invariant_drop
	print "Number of loci with GC content within the threshold:", counter_GChigh #FFM
	print "Number of loci with 3 alleles in at least 1 individual:", counter_3all_drop
	print "Number of loci deleveraged in Stacks:", counter_deleveraged
	print "Number of loci showing observed heterozygosity higher than allowed:", counter_het_dropped
	print "Number of loci with Ns in any of the samples:", counter_Ns_dropped 
	print "Number of usable loci:", len(cat_id_list)
	print "Average sample size per locus before:", sum(list_ind_by_locus)/len(list_ind_by_locus),', and after filtering:', sum(list_sample_size)/len(list_sample_size)


	return cat_id_list, dropped_loci, global_position_list, GC_content_list, list_ind_by_locus, list_sample_size


# Aqui va el script FFM

import sys
import numpy
from decimal import *
import matplotlib.pyplot as plt
from pylab import subplot
from pylab import subplots_adjust
import argparse as ap

parser = ap.ArgumentParser()

parser.add_argument('-i', '--infile', help='Provide a tsv file as that created by export_sql.pl in the Stacks (Catchen et al 2013) package', required=True, type=str)
# FFM
parser.add_argument('-p1', '--position_min', help='Minimum position in the sequence to retain a SNPs (SNPs before this position will be discarded)', required=True, type=int, default=0)
# FFM
parser.add_argument('-p2', '--position_max', help='Maximum position in the sequence to retain a SNPs (SNPs after this position will be discarded)', required=True, type=int)
parser.add_argument('-t', '--heterozygosity', help='Maximum locus heterozygosity to keep a locus in the database  (between 0 and 1)', required=True, type=float)
parser.add_argument('-c1', '--GC_content_min', help='Minimum locus GC content to keep a locus in the database (between 0 and 1)', required=True, type=float, default=0)
parser.add_argument('-c2', '--GC_content_max', help='Maximum locus GC content to keep a locus in the database (between 0 and 1)', required=True, type=float)
parser.add_argument('-m', '--missing', help='Minimum number of individuals per locus', required=True, type=float)
args = parser.parse_args()

#launch the function
snp_pos_min = args.position_min # FFM
snp_pos_max = args.position_max #FFM
het_max_allowed = args.heterozygosity
GC_cont_min = args.GC_content_min
GC_cont_max = args.GC_content_max
input_file = args.infile 
max_missing_ind = args.missing

cat_id_list, dropped_loci, global_position_list, GC_content_list, list_ind_by_locus, list_sample_size = loci_selector(input_file,snp_pos_min,snp_pos_max,het_max_allowed,GC_cont_min,GC_cont_max,max_missing_ind)


#write the catalog ID of the selected loci to a file ready to be used in grep
output_file = "ID_list_"
output_file += input_file
out = open(output_file, "w")
for IDs in cat_id_list:
	out.write('^'+IDs+'\s'+'\n')
out.close()

#write the catalog ID of the dropped loci to a file ready to be used in grep
output_file = "dropped_loci_"
output_file += input_file
out = open(output_file, "w")
for IDs in dropped_loci:
	out.write('^'+IDs+'\s'+'\n')
out.close()

#histogram the SNPs position across all retained loci, GC content and sample size distribution
print 'Average position of the SNPs along the reads is being plotted as well as GC content across selected loci'

subplot(4,1,1)
plt.hist(numpy.asarray(global_position_list, dtype='int'), bins=range(0,snp_pos_max+1), normed=False)
plt.title('SNPs position across loci in '+input_file)
plt.xlabel('position')

subplot(4,1,2)
plt.hist(numpy.asarray(GC_content_list, dtype='float'), bins=20, normed=False)
plt.title('Mean GC content across loci')
plt.xlim(0.2,0.8)
plt.xlabel('GC')

subplot(4,1,3)
plt.hist(numpy.asarray(list_ind_by_locus, dtype='int'), bins=range(min(list_ind_by_locus)-1,max(list_ind_by_locus)+1), normed=False)
plt.title('Average number of individuals across loci before filtering')
plt.xlabel('Individuals')

subplot(4,1,4)
plt.hist(numpy.asarray(list_sample_size, dtype='int'), bins=len(set(list_sample_size)), normed=False)
plt.title('Average number of individuals across loci after filtering')
plt.xlabel('Individuals')

subplots_adjust(hspace=.4)

plt.show()



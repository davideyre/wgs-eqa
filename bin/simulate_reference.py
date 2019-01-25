#!/usr/bin/env python3


## generate a simulated reference sequence with mutations, recombination and indels
## David Eyre, david.eyre@bdi.ox.ac.uk
## 25 January 2019

## Could potentially use to map reads back
## If wish to only recover SNPs then need to output a list of variants

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import sys
import numpy as np
import argparse
import uuid


def simulate_ref (ref_file, output, mutation_n, subs_from_recomb, recomb_length_mean, prop_of_sites_diff_in_recomb, deletion_n, deletion_length_mean, insertion_n, insertion_length_mean):

	#get the expected number of recombination events
	recomb_events = round(subs_from_recomb / (prop_of_sites_diff_in_recomb * recomb_length_mean))
	
	# read in reference and length over all chromosomes
	ref = [chr for chr in SeqIO.parse(ref_file, 'fasta')]
	ref_len = sum(len(chr) for chr in ref)
	
	#base frequencies
	concat_seq = "".join(chr.seq._data for chr in ref)
	bc = Counter(concat_seq)
	bases = ['A', 'C', 'G', 'T']
	pctA = round(float(bc['A'])/ref_len,5)
	pctC = round(float(bc['C'])/ref_len,5)
	pctG = round(float(bc['G'])/ref_len,5)
	pctT = 1.00 - pctA - pctC - pctG
	bases_p = [pctA, pctC, pctG, pctT]
	
	#get event rates
	mutation_per_site = float(mutation_n) / ref_len
	recomb_events_per_site = float(recomb_events) / ref_len
	deletion_events_per_site = float(deletion_n) / ref_len
	insertion_events_per_site = float(insertion_n) / ref_len
	
	#events
	event_per_site = mutation_per_site + recomb_events_per_site + deletion_events_per_site + insertion_events_per_site
	event_list = ['mutation', 'recomb', 'deletion', 'insertion']
	event_list_p = [e/event_per_site for e in [mutation_per_site, recomb_events_per_site, deletion_events_per_site, insertion_events_per_site]]
	
	#output new_ref
	new_ref = []
	mutation_counter = 0
	recomb_counter = 0
	recomb_subs = 0
	del_counter = 0
	del_bases = 0
	insert_counter = 0
	insert_bases = 0
	
	#ouput list of variants
	output_list_file = '%s.txt'%output
	o = open(output_list_file, 'w')
	o.write('variant_type\tchromosome\tsite\tnew_site\toriginal\tnew\n')
	
	for chr in ref:
		#for each chromosome / plasmid in reference
		seq_data = [b for b in chr.seq]
		chr_len = len(chr)
		insert_list = []
		insert_bases_chr = 0
		site_offset = 0
		
		#get first site with an event
		site = np.random.geometric(event_per_site)
		
		while (site < chr_len):
			
			#select an event type
			event_type = np.random.choice(event_list, 1, p=event_list_p)[0]
			if event_type=='mutation':
				#a mutation has occurred
				mutation_counter += 1
				
				#simulate the new base
				new_base = np.random.choice(bases, 1, p=bases_p)[0]
				o.write('mutation\t%s\t%s\t%s\t%s\t%s\n'%(chr.id, site, site+site_offset, seq_data[site], new_base))			
				#add mutation to new reference
				seq_data[site] = new_base
				
			
			elif event_type == "deletion":
				#a deletion event has occurred
				del_counter += 1
				#get number of bases affected
				del_len = np.random.geometric(1/deletion_length_mean)
				#check not beyond end of chr
				if site+del_len > chr_len:
					#if it is set the end of the recombination tract to the end of the chr
					end = chr_len
				else:
					end = site+del_len
				#set deletion
				for b in range(site, end):
					o.write('deletition\t%s\t%s\t%s\t%s\t-\n'%(chr.id, b, b+site_offset, seq_data[b]))
					site_offset -= 1	
					seq_data[b] = '-'
				del_bases += end-site
			
			elif event_type == "insertion":
				#insertion event
				insert_len = np.random.geometric(1/insertion_length_mean)
				insert = "".join(np.random.choice(bases, insert_len, p=bases_p))
				insert_list.append((site+insert_bases_chr, insert))
				o.write('insertion\t%s\t%s\t%s\t%s\t%s\n'%(chr.id, site, site+site_offset, seq_data[site], seq_data[site] + insert))
				insert_counter += 1

				insert_bases_chr += insert_len
				site_offset += insert_len
				
			else:
				#recombination has occurred
				recomb_counter += 1
				
				#simulate a recombination length
				event_len = np.random.geometric(1/recomb_length_mean)
				#check event end isn't beyond end of chr
				if site+event_len > chr_len:
					#if it is set the end of the recombination tract to the end of the chr
					end = chr_len
				else:
					end = site+event_len
				
				#set the mutations in the recombination tract
				for b in range(site, end):
					if np.random.random() < prop_of_sites_diff_in_recomb:
						new_base = np.random.choice(bases, 1, p=bases_p)[0]
						o.write('recombination\t%s\t%s\t%s\t%s\t%s\n'%(chr.id, b, b+site_offset, seq_data[b], new_base))	
						seq_data[b] = new_base
						recomb_subs += 1
			
			#simulate number of sites to next event
			site = site + np.random.geometric(event_per_site)
		
		#once complete write the seq_data back
		
		#insert insertions - the insertion follows the base identified by siite
		insert_bases += insert_bases_chr #keep track of the overall number of inserted bases
		last_i = 0
		seq_data_i = []
		for i in insert_list:
			insert_seq = [b for b in i[1]]
			if (last_i==0):
				seq_data_i = seq_data[0:i[0]+1] + insert_seq
			else:
				seq_data_i += seq_data[last_i+1:i[0]+1] + insert_seq
			last_i = i[0]
		
		seq_data_i += seq_data[last_i+1:]
		
		seq_data_joined = "".join([b for b in seq_data_i if b!="-"])
		desc_str = "| m=%s sr=%s rl=%s p=%0.3f d=%s dl=%s i=%s il=%s"%(
							mutation_n, subs_from_recomb, recomb_length_mean, prop_of_sites_diff_in_recomb,
							deletion_n, deletion_length_mean, insertion_n, insertion_length_mean)
		seq_record = SeqRecord(seq=Seq(seq_data_joined), id=chr.id, 
								name=chr.name, description=chr.description + desc_str)
		new_ref.append(seq_record)
	
	#write the mutated reference file
	out_fa = '%s.fa'%output
	SeqIO.write(new_ref, out_fa, 'fasta')
	new_ref_len = sum(len(chr) for chr in new_ref)
	o.close()
	
	#report number of events of each type
	sys.stdout.write('\nSimulated reference written to %s\n'%out_fa)
	sys.stdout.write('List of variants written to %s\n'%output_list_file)
	sys.stdout.write('There are %s mutations\n'%mutation_counter)
	sys.stdout.write('There are %s recombinations resulting in %s substitutions\n'%(recomb_counter, recomb_subs))
	sys.stdout.write('There are %s deletions involving a total of %s bases\n'%(del_counter, del_bases))
	sys.stdout.write('There are %s insertions involving a total of %s bases\n'%(insert_counter, insert_bases))
	sys.stdout.write('Original length: %s\tNew length %s\t Difference %s\n'%(ref_len, new_ref_len, new_ref_len-ref_len))
	sys.stdout.write('Done.\n\n')

if __name__ == '__main__':
	
	#parse inputs
	parser = argparse.ArgumentParser(description='Simulate a new reference from an existiing reference with mutation, recombination, indels')

	parser.add_argument('-r', '--refile', dest='ref_file', action='store', help='path to existing reference file', required=True)
	parser.add_argument('-o', '--output', dest='output', action='store', help='path stem for new reference file', required=True)
	
	parser.add_argument('-m', '--mutation_n', dest='mutation_n', action='store', type=int, default=500,
						help='expected number of mutation events (including same base to same base)')
	parser.add_argument('-s', '--subs_from_recomb', dest='subs_from_recomb', action='store', type=int, default=300,
						help='expected number of substitutions from recombination events (including same base to same base)')
	parser.add_argument('-l', '--recomb_length_mean', dest='recomb_length_mean', action='store', type=int, default=1000,
						help='mean recombination fragment length')
	parser.add_argument('-p', '--prop_of_sites_diff_in_recomb', dest='prop_of_sites_diff_in_recomb', action='store', type=float, default=0.04,
						help='proportion of sites substituted in recombination fragment (includng same to same)')
	parser.add_argument('-d', '--deletion_n', dest='deletion_n', action='store', type=int, default=100,
						help='mean number of deletion events')
	parser.add_argument('-e', '--deletion_length_mean', dest='deletion_length_mean', action='store', type=float, default=1.5,
						help='mean length of deletion event (geometrically distributed)')
	parser.add_argument('-i', '--insertion_n', dest='insertion_n', action='store', type=int, default=100,
						help='mean number of insertion events')
	parser.add_argument('-n', '--insertion_length_mean', dest='insertion_length_mean', action='store', type=float, default=1.5,
						help='mean length of insertion event (geometrically distributed)')
	args = parser.parse_args()
	
	ref_file = args.ref_file
	output = '%s_%s'%(args.output, str(uuid.uuid4())[0:8])
	mutation_n = args.mutation_n
	subs_from_recomb = args.subs_from_recomb
	recomb_length_mean = args.recomb_length_mean
	prop_of_sites_diff_in_recomb = args.prop_of_sites_diff_in_recomb
	deletion_n = args.deletion_n
	deletion_length_mean = args.deletion_length_mean
	insertion_n = args.insertion_n
	insertion_length_mean = args.insertion_length_mean
	
	#run the simulation
	simulate_ref (ref_file, output, mutation_n, subs_from_recomb, recomb_length_mean, prop_of_sites_diff_in_recomb, deletion_n, deletion_length_mean, insertion_n, insertion_length_mean)
	
	#example call
	#./simulate_reference.py -r ref/R00000003.fasta -o example_output/R00000003 -m 500 -s 300 -l 1000 -p 0.04 -d 100 -e 1.5 -i 100 -n 1.5
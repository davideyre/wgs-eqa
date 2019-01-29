#!/usr/bin/env python3

from Bio import SeqIO
import uuid, argparse, os, sys, json

def parse_infile(infile):
	#read in CSV file
	# columns - stem to simulated .fa/.json/.txt files, concensus fasta file from mapping
	filelist = []
	with open(infile) as o:
		o.readline() #skip header
		for l in o:
			l = l.strip().split(',')
			filelist.append(l)
	return filelist


def compare_sim(refFile, simList, outFile):
	#read in reference file
	ref = [s for s in SeqIO.parse(refFile, 'fasta')]
	refLen = sum([len(s) for s in ref])
	
	#open output file
	outHandle = open(outFile, 'w')
	
	#iterate all over simulations
	for sim in simList:
		testFile = sim[1]
		test = [s for s in SeqIO.parse(testFile, 'fasta')]
		
		#read in the list of true variants
		truthid = os.path.basename(sim[0])
		truthFile = '%s.txt'%sim[0]
		truth = {}
		var_count = {'mutation': 0, 'recombination': 0}
		with open(truthFile, 'r') as o:
			header = o.readline() #skip the header
			for l in o:
				l = l.strip().split()
				if l[4] != l[5]: #ignore same base to same base changes
					if l[3] in ('mutation', 'recombination'): #don't track inserts and deletions
						if l[0] not in truth.keys():
							truth[l[0]] = {}
						truth[l[0]][l[1]] = {'type': l[3], 'original': l[4], 'variant': l[5], 'found': False}
						var_count[l[3]] += 1 #record the total number of variants to find
		
		#read in the metadata about the simulation
		jsonFile = '%s.json'%sim[0]
		with open(jsonFile) as o:
			truthDict = json.load(o)
		
		tp_count = {'mutation': 0, 'recombination': 0}
		fp_count = 0
		# check which variants found
		for i, chr in enumerate(ref):
			for site, (b_ref, b_test) in enumerate(zip(ref[i], test[i])):
				if b_ref != b_test and b_test in 'ACGT':
					#check for false positive
					if str(site) not in truth[chr.id].keys():
						#print ('FALSE POSITIVE: %s'%site)
						fp_count += 1
						#print(truthid, testid, chr.id, site, b_ref, b_test)
					#check for match
					else:
						if b_test == truth[chr.id][str(site)]['variant']:
							#print ('TRUE POSITIVE: %s (%s)'%(site , truth[chr.id][str(site)]['type']))
							truth[chr.id][str(site)]['found'] = True
							tp_count[truth[chr.id][str(site)]['type']] += 1
		
		var_total = 0
		tp_total = 0
		for var_type in var_count.keys():
			var_total += var_count[var_type]
			tp_total += tp_count[var_type]
		
		sens = (100 * float(tp_total)) / var_total #percentage sensitivity
		spec = float(fp_count) / refLen * 1e6 #errors per 1 Mb of reference
		
		outStr = '%s\t%0.2f\t%0.2f\t%s\n'%(truthid, sens, spec, truthDict['parms'])
		sys.stdout.write(outStr)
		outHandle.write(outStr)
	outHandle.close()

if __name__ == '__main__':
	
	#parse inputs
	parser = argparse.ArgumentParser(description='Compare pipeline performance with simulation')
	parser.add_argument('-i', '--in_file', dest='inFile', action='store', 
						help='csv file containing details of input files', required=True)
	parser.add_argument('-r', '--ref_file', dest='refFile', action='store', 
						help='reference fasta file used for mapping', required=True)						
	parser.add_argument('-o', '--output_file', dest='outFile', action='store', 
						help='path for summary file', required=True)
	args = parser.parse_args()
	
	#parse list of simulations
	simList = parse_infile(args.inFile)
	
	#compare the simulations to the original reference
	compare_sim(args.refFile, simList, args.outFile)
	
	#example call - 
	# bin/compare_simulations.py -r ref/R00000003.fasta -i ../pipeline-test/sim-compare.csv -o ../pipeline-test/sim_compare_results.txt

			
	
	
		

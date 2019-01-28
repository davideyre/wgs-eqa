#!/usr/bin/env python

# collect simulations and write an output file for use with bug-flow

import uuid
import argparse
import os

if __name__ == '__main__':
	
	#parse inputs
	parser = argparse.ArgumentParser(description='Collect simulations into CSV')
	parser.add_argument('-p', '--path', dest='path', action='store', 
						help='path to simulated reference', required=True)
	parser.add_argument('-o', '--output_file', dest='out_file', action='store', 
						help='path for summary file', required=True)
	args = parser.parse_args()
	
	
	#list all the files in the input director
	dir_list = os.listdir(args.path)
	
	with open(args.out_file, 'w') as o:
		o.write('sampleid,uuid,fq1,fq2\n')	
		for f in dir_list:
			if f.endswith(".fa"):
				id = f.split('.')[0]
				fq1 = 'r%s.1.fq.gz'%id
				fq2 = 'r%s.2.fq.gz'%id
				uuid = str(uuid.uuid4())
				o.write('%s,%s,%s,%s\n'%(id, uuid, fq1, fq2))
#!/usr/bin/env nextflow

// parameters 
params.outputPath = "example_output"
params.refFile = "ref/R00000003.fasta"
params.iterations = 10

// initial logging
log.info "\n" 
log.info "WGS-EQA simulate -- version 0.1"
log.info "Input reference file   :  ${params.refFile}"
log.info "Output path            :  ${params.outputPath}"
log.info "\n"

// rename input parameters
refFasta = file(params.refFile)
outputPath = file(params.outputPath)
iterations = params.iterations

// set up initial channel with sets of 4 values - mutations, substitutions from recombinations, deletions, insertions
Channel
    .from( [500, 0, 0, 0], [500, 500, 0, 0], [500, 500, 200, 200], [2500, 2500, 500, 500], [25000, 25000, 1000, 1000] )
    .set { simSettings }


// make simulated references
process simReference {
  
    input:
        file refFasta
        file outputPath
        set (m, s, d, i) from simSettings
        each iter from 1..iterations
	
	output:
		file "*"
		file "*.fa" into faList
	
	tag {refFasta}
	publishDir "$outputPath", mode: 'copy'

	//for now fix the length of the recombination fragment, the diversity of the 
	//recombination fragement, and the deletion and insertions size distribution:
	// -l 500 -p 0.05  -e 1.5 -n 1.5
    """
    simulate_reference.py -r $refFasta \
    	-o . \
    	-m $m -s $s -d $d -i $i \
    	-l 500 -p 0.05  -e 1.5 -n 1.5
    """
}

//simulate reads for all 

process simReads {

	input:
		file fasta from faList
	
	output:
		file "*.gz"
	
	publishDir "$outputPath", mode: 'copy'
	
	//genome length = 4.3e6 * 50 (depth) / 150 (read length) = 1433333
	//-e 0.001 expected sequencing error rate
	// -1 150 -2 150 150bp reads
	// -r 0 -R 0 don't simulate more mutations or indels
	// - S random seed = 1
	
	//prefix output files with r for rSR - reads from Simulated Reference
	"""
	wgsim -N 1433333 -e 0.001 -1 150 -2 150 -r 0 -R 0 -S 1  ${fasta} r${fasta.baseName}.1.fq r${fasta.baseName}.2.fq
	gzip *.fq
	"""
	
}
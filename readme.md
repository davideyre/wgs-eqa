# WGS EQA
A set of tools for external quality assurance for pathogen sequencing

## Dependencies
 - Python 3
 - Numpy
 - Biopython

## Contents
1. Simulate

## 1.Simulate
Generates a simulated reference sequence from an existing reference sequence. Allows for mutations, recombination events, and indels to be simulated.

Example usage
```
bin/simulate_reference.py -r ref/R00000003.fasta -o example_output/R00000003 -m 500 -s 300 -l 1000 -p 0.04 -d 100 -e 1.5 -i 100 -n 1.5
```

Example usage, simulating mutations only
```
bin/simulate_reference.py -r ref/R00000003.fasta -o example_output/R00000003 -m 500 -s 0 -d 0  -i 0
```

Obtain help
```
bin/simulate_reference.py -h
usage: simulate_reference.py [-h] -r REF_FILE -o OUTPUT [-m MUTATION_N]
                             [-s SUBS_FROM_RECOMB] [-l RECOMB_LENGTH_MEAN]
                             [-p PROP_OF_SITES_DIFF_IN_RECOMB] [-d DELETION_N]
                             [-e DELETION_LENGTH_MEAN] [-i INSERTION_N]
                             [-n INSERTION_LENGTH_MEAN]

Simulate a new reference from an existiing reference with mutation,
recombination, indels

optional arguments:
  -h, --help            show this help message and exit
  -r REF_FILE, --refile REF_FILE
                        path to existing reference file
  -o OUTPUT, --output OUTPUT
                        path stem for new reference file
  -m MUTATION_N, --mutation_n MUTATION_N
                        expected number of mutation events (including same
                        base to same base)
  -s SUBS_FROM_RECOMB, --subs_from_recomb SUBS_FROM_RECOMB
                        expected number of substitutions from recombination
                        events (including same base to same base)
  -l RECOMB_LENGTH_MEAN, --recomb_length_mean RECOMB_LENGTH_MEAN
                        mean recombination fragment length
  -p PROP_OF_SITES_DIFF_IN_RECOMB, --prop_of_sites_diff_in_recomb PROP_OF_SITES_DIFF_IN_RECOMB
                        proportion of sites substituted in recombination
                        fragment (includng same to same)
  -d DELETION_N, --deletion_n DELETION_N
                        mean number of deletion events
  -e DELETION_LENGTH_MEAN, --deletion_length_mean DELETION_LENGTH_MEAN
                        mean length of deletion event (geometrically
                        distributed)
  -i INSERTION_N, --insertion_n INSERTION_N
                        mean number of insertion events
  -n INSERTION_LENGTH_MEAN, --insertion_length_mean INSERTION_LENGTH_MEAN
                        mean length of insertion event (geometrically
                        distributed)
```

Having run the python script to simulate a new reference, wgsim can be used to simulate read data. wgsim can be downloaded and compiled from source:
```
git clone https://github.com/lh3/wgsim.git
cd wgsim
gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm
```

Or installed via bioconda:

```
conda config --add channels bioconda
conda install wgsim
```


An example of this approach, as a NextFlow workflow can be found in simulate.nf. 

An alternative approach is to map actual sequence reads from sequencing the original reference back to the simulated reference with variants.


David Eyre
david.eyre@bdi.ox.ac.ul
25 January 2019
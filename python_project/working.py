#!/usr/bin/env python

"""
This is a working version of the script find_best_genomebins.py

The script takes a user-input assembly file and a user-input mapping file generated from 
binning (of contigs to bins (must be an integer)) 
that can be either tab or comma seperated.

The script can be used with output from any binning method as long as you have the original 
metagenome assembly used and a file mapping the contigs 
to genome bin numbers.

Some examples of possible binning methods (either user or software generated bins) include 
CONCOCT, PhymmBL and ESOM (user-selected bins). This script is modellingafter using CONCOCT 
output data but it should be possible to use output from any binning method if it is in the 
specified format.

The script will output descriptive data about the different input bins both graphically 
(as a histogram distribution)in a pdf file and as a csv table.

Descriptive data includes: 
- GC content
- genome bin length (bp)
- N50 statistic (bp)
- *completness (in future expansion)
- *redundancy (in future expansion)

In a given output folder (which if it does not exist will be created) the script will also 
generate sub-folders for each of the genome bins containing a file with contigs assigned to 
each in fasta format (file.fna).  

*FUTURE EXPANSION
- use checkM to get completness and redundancy
- the user can then specify thresholds for both variables (although there are default values)
- if the user does not specify the flag "--all" only bins which pass the given threshold will
be output and included in the descriptibe data

HELP


"""

#To use this script matplotlib and Biopython need to be installed
import sys
import os
import re
#import argparse
#import matplotlib.pyplot as plt
from Bio import SeqIO
#from Bio.SeqUtils import GC
#from Bio.Alphabet import IUPAC

#Biopython methods to use
#Bio.SeqIO.to_dict()

#############################################################################################

class Assembly(object):
	"""
	This object contains the attributes of the path to the input assembly fasta file and the 
	input mapping file of contig assignment to bins (output by hand or from binning software)

	It includes a method extract_bins for extracting the sequences of contigs associated
	with each bin.
	"""
	def __init__(self, fasta_file, binning_file):

		if not os.path.exists(fasta_file):
			raise IOError("This file does not exist")

		handle = open(fasta_file,'rU')
		record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
		self.fasta_file = record_dict

		if not os.path.exists(binning_file):
			raise IOError("This file does not exist")

		handle = open(binning_file,'rU')
		parsed_bin_file = {}
		for i in handle:
			contig, bin_number = i.strip('\n').split(',')
			bin_number = int(bin_number)
			if bin_number in parsed_bin_file.keys():
				parsed_bin_file[bin_number].append(contig)
			else:
				parsed_bin_file[bin_number] = []
				parsed_bin_file[bin_number].append(contig)
		self.binning_file = parsed_bin_file
	
	def extract_bins(self):
		pass
		


trial = Assembly("test.fasta", "mapping.txt")
print trial.fasta_file
print trial.binning_file




#class Bin(object):
	#def __init__(self, bin_number, contigs):
	#	self.bin_number = bin_number
	#	self.contigs = contigs
	#	self.num_contigs = num_contigs
	#	self.GC_content_perc = GC_content_perc
	#	self.bin_length_bp = bin_length_bp
	#	self.N50 = N50
		#self.completness = completness
		#self.redundancy = redundancy
		#self.best_bin = best_bin

#class Contig(object):
	#def __init__(self):
	#	pass

#class Checkm(Bin):
	#This will be added at a later point to calculate redundancy and completness of each bin
#	def __init__(self):
#		pass

#############################################################################################

#parser = argparse.ArgumentParser(prog='find_best_genomebins.py', 
#	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
#	description='Generate descriptive information and fasta files for genome bins from an \
#	assembly fasta file and associated binning output (contigs mapped to bin numbers)')

#parser.add_argument('-assem', '--assembly_path', required=True, help = 'Path to the metagenome assembly file in \
#	fasta format with the extension .fna or .fasta')

#parser.add_argument('-map', '--contig_to_bin_map_path', required=True, help = 'Path to the mapping file containing contigs with \
#	associated bin assignments in the format "contig_name, bin_#", with the "bin_#" as an \
#	integer. The contig and bin can be either comma or tab seperated. It should be a text \
#	file with either a .txt extension or no extension')

#parser.add_argument('-o', '--output_folder', required=True, help = 'Completness threshold a genome bin has to pass to be \
#	output and described. As a proportion (float from and including 0 to 1)')

#parser.add_argument('--compl', type=float, default=0.7, help = 'Completness threshold a genome bin has to pass to be \
#	output and described. As a proportion (float from and including 0 to 1)')

#parser.add_argument('--redun', type=float, default=1.2, help = 'Redundancy threshold a genome bin has to pass to be \
#	output and described. As a proportion (float >= 1)')

#parser.add_argument('--all', action='store_true', help = 'Include this flag if you want all bins output and described \
#	instead of only those that pass the given completness and redundancy thresholds (the default)')

#parser.add_argument('--t', type=int, default=1, help = 'Number of threads for checkM to use when calculated the \
#	completness and redundancy of each bin')
#args = parser.parse_args()

#print "The path given for the assembly: " + args.assembly_path
#print "The path given for the mapping file: " + args.map_path
#print "The completness threshold is " + args.compl
#print "The redundancy threshold is: " + args.redun
#print "All bins will be output: " + args.all
#print "Number of threads used by CheckM: " + args.t
#print "Files will be output in this directory: " + args.o
#parser.print_help()

#############################################################################################

#Acutal implementation
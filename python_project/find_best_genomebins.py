#!/usr/bin/env python

"""
This is a working version of the script find_best_genomebins.py

The script takes a user-input assembly file and a user-input mapping 
file generated from binning (of contigs to bins) 
that needs to be comma seperated.

The script can be used with output from any binning method as long as 
you have the original metagenome assembly used and a file mapping the 
contigs to genome bin numbers.

Some examples of possible binning methods (either user or software 
generated bins) include CONCOCT, PhymmBL and ESOM (user-selected bins). 
This script is modellingafter using CONCOCT output data but it should 
be possible to use output from any binning method if it is in the 
specified format.

The script will output descriptive data about the different input bins 
both visually (as an interactive plot saved in a html file) and as a 
csv table.

Descriptive data includes: 
- Number of contigs
- N50 statistic (bp)
- GC content
- genome bin length (bp)
- *completness (in future expansion)
- *redundancy (in future expansion)

In a given output folder (which if it does not exist will be created) 
the script will also generate sub-folders for each of the genome bins 
containing a file with contigs assigned to each in fasta format 
(file.fasta). 

*FUTURE EXPANSION
- A log file will be created
- use checkM to get completness and redundancy
- the user can then specify thresholds for both variables 
(although there are default values)
- if the user does not specify the flag "--all" only bins which pass 
the given threshold will be output and included in the descriptibe data

HELP:

usage: find_best_genomebins.py [-h] -assem ASSEMBLY_PATH -map CONTIG_TO_BIN_MAP_PATH -o
                  OUTPUT_FOLDER

Generate descriptive information and fasta files for genome bins from an
assembly fasta file and associated binning output (contigs mapped to bin
numbers)

optional arguments:
  -h, --help            show this help message and exit
  -assem ASSEMBLY_PATH, --assembly_path ASSEMBLY_PATH
                        Path to the metagenome assembly file in fasta format
                        with the extension .fna or .fasta (default: None)
  -map CONTIG_TO_BIN_MAP_PATH, --contig_to_bin_map_path CONTIG_TO_BIN_MAP_PATH
                        Path to the mapping file containing contigs with
                        associated bin assignments in the format "contig_name,
                        bin_#", with the "bin_#" as an integer. The contig and
                        bin can be either comma or tab seperated. It should be
                        a text file with either a .txt extension or no
                        extension (default: None)
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Path to the output folder, will be created if it does
                        not already exist (default: None)
"""

##############
#To use this script Biopython and bokeh need to be installed (easiest using pip)

#General modules
import sys
import os
import re
import logging
import argparse
import numpy as np
import pandas as pd

#Modules for sequence analysis
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
from Bio.Alphabet import IUPAC

#Modules for data visualization
from bokeh.charts import Scatter, save
from bokeh.charts import output_file as out_file

#############################################################################################


class Assembly(object):


	"""
	This class contains the attributes of the path to the input assembly
	fasta file, the input mapping file of contig assignment to bins and
	the output folder path to which new files and folders will be written.

	The parse_fasta_file and parse_binning_file methods parse the data 
	in the two input files into dictionary structures.

	The make_Bin_instances method makes new instances of the second 
	class Bin.

	The methods table_bin_characteristics and plot_bin_characteristics 
	output descriptive data in a csv and html file, respectively, about
	 the different bins from the assembly.

	The extract_bins method creates folders for each bin inside the 
	specified output folder and writes a fasta file containing contigs 
	that were mapped to each bin inside the bin folder.
	"""


	def __init__(self, fasta_file, binning_file, output_folder):

		"""
		These attributes must be passes to a new instance:

		fasta_file = the path to the metagenomic assembly fasta file
		binning_file = the path to the mapping file of contigs to bins
		output_folder = path to a directory for bins to be output
		(if not already existing it will be created)

		Attributes created using methods:

		bins_dict = a dictinary of bins, containing a dictionary of 
		contigs as SeqRecord objects (containing sequence information) 
		bins = a list of Bin class instances
		bin_characteristics = a table containing characteristics 
		(calculated in the Bin class) of each Bin instance
		"""

		if not os.path.exists(fasta_file):
			raise IOError("This file does not exist")
		self.fasta_file = fasta_file

		if not os.path.exists(binning_file):
			raise IOError("This file does not exist")
		self.binning_file = binning_file
		
		if not os.path.exists(output_folder):
   				os.makedirs(output_folder)
   		self.output_folder = output_folder

		self.bins_dict = {}
		self.bins = []
		self.bin_characteristics = None

	def parse_fasta_file(self):

		"""
		The fasta file is parsed into biopython SeqRecord object in 
		a SEQIO dictionary

		Keys are contig names
		and values are the SeqRecord object
		use .seq to get the nuceltoide sequence
		"""

		handle = open(self.fasta_file,'rU')
		record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta", IUPAC.unambiguous_dna))
		self.fasta_file = record_dict
		handle.close()

	def parse_binning_file(self):

		"""
		The binning file is parsed into a dictionary of bin
		numbers, each containing a list of contig names.

		File must be in comma-seperated format 
		*Later Expansion:
		will be expanded to allow tab-seperated format at a later point
		"""

		handle = open(self.binning_file,'rU')
		parsed_bin_file = {}
		for i in handle:
			contig, bin_number = i.strip('\n').split(',')
			if bin_number in parsed_bin_file.keys():
				parsed_bin_file[bin_number].append(contig)
			else:
				parsed_bin_file[bin_number] = []
				parsed_bin_file[bin_number].append(contig)
		self.binning_file = parsed_bin_file
		handle.close()

	def make_Bin_instances(self):

		"""
		This method creates instances of the class Bin.

		It first creates a dictionary (self.bins_dict) with each bin as
		a key and values as an imbedded dictionary with contigs as keys
		for SeqRecord objects containing associated sequences.

		Then a list (self.bins) of Bin class instances with each bin as
		the attribute "bin_number" and associated contigs as the 
		attribute "contigs" is created.
		"""

		for genome_bin in self.binning_file:
			self.bins_dict[genome_bin] = {}
			for contig in self.binning_file[genome_bin]:
				if contig == self.fasta_file[contig].id:
					self.bins_dict[genome_bin][contig] = self.fasta_file[contig]
		for genome_bin in self.bins_dict:
			self.bins.append(Bin(genome_bin, self.bins_dict[genome_bin]))

	def table_bin_characteristics(self):

		"""
		A series of lists of bin characteristics corresponding to Bin 
		class attributes (bin_name, number of contigs, N50 statistic, 
		GC content and genome bin length) are created. 

		These are combined into a pandas dataframe 
		(self.bin_characteristics) using bin numbers as an index.

		A tab-seperated csv file is created from the dataframe.
		"""

		index_number = [genome_bin.bin_number for genome_bin in self.bins]
		bin_number_list = [('bin_' + genome_bin.bin_number) for genome_bin in self.bins]
		number_contigs_list = [genome_bin.number_contigs for genome_bin in self.bins]
		N50_list = [genome_bin.N50 for genome_bin in self.bins]
		GC_list = [genome_bin.GC for genome_bin in self.bins]
		bin_length_bp_list = [genome_bin.bin_length_bp for genome_bin in self.bins]
		
		d = {'bin_name' : pd.Series(bin_number_list, index=index_number),
		'number_of_contigs' : pd.Series(number_contigs_list, index=index_number),
 		'N50_statistic' : pd.Series(N50_list, index=index_number),
 		'GC_content' : pd.Series(GC_list, index=index_number),
 		'length_bp' : pd.Series(bin_length_bp_list, index=index_number)
 		}

		df = pd.DataFrame(d)
		cols = ['bin_name', 'number_of_contigs', 'N50_statistic', 'GC_content', 'length_bp']
		self.bin_characteristics = df[cols]
		self.bin_characteristics.to_csv((self.output_folder + "/" + "bin_characteristics.csv"),
										 sep='\t')

	def plot_bin_characteristics(self):
		
		"""
		GC content and genome bin length are plotted in a scatterplot
		from the self.bin_characteristics attribute using bokeh.

		The plot is output as an html file and is interactive.
		By however over any datapoint (each representing a bin),
		other descriptive information about the bin and it's name
		can be seen.
		"""

		tooltips=[
    		('Bin Name', '@bin_name'),
    		('Number of Contigs', '@number_of_contigs'),
    		('N50 Statistic', '@N50_statistic'),
    		('GC Content (%)', '@GC_content'),
    		('Length (bp)', '@length_bp')
    		]

		p = Scatter(self.bin_characteristics, 
			x='length_bp', 
			y='GC_content', 
			title="GC Content vs Length of Genome Bins",
            xlabel="GC Content (%)", 
            ylabel="Length (bp)",
            tooltips=tooltips
            )

		destination = self.output_folder + "/" + "bin_characteristics.html"
		out_file (destination)
		save(p)
	
   	def extract_bins(self): 

   		"""
		This method uses the parsed metagenome assembly and mapping 
		file to output folders for every bin, each containing a fasta
		file with contigs mapped to it.

		*Later Expansion:
		Will be altered later to only extract bins where best_bin is 
		True (based on completness and redundancy threshold, unless 
		the user has input the -all flag to have all sequences output)
		"""

		for genome_bin in self.binning_file:
			if not os.path.exists(self.output_folder + "/" + ("bin_%s" % (genome_bin))):
				os.makedirs(self.output_folder + "/" + ("bin_%s" % (genome_bin)))
		for genome_bin in self.binning_file.keys():
			output_file = self.output_folder + "/" + ("bin_%s" % (genome_bin)) + "/" + ("bin_%s" % (genome_bin)) + ".fasta"
			contig_names = self.binning_file[genome_bin]
			record = []
			for name in contig_names:
				record.append(self.fasta_file[name])
			SeqIO.write(record, output_file, "fasta")
			print "%s folder made and contigs extracted" % (genome_bin)
		print "Bin extraction complete"


class Bin(object):


	"""
	This class contains attributes corresponding to characteristics of 
	each bin as well as bin numbers and SeqRecord objects (contigs) 
	associated with each bin. 

	The calculate_bin_descriptors method calculates basic descriptive 
	information about each bin given bin_numbers and contigs.
	"""


	def __init__(self, bin_number, contigs):

		"""
		These attributes must be passes to a new instance:

		bin_number = the number of the bin, an integer as a string
		contigs = a dictionary including contigs mapped to the bin
		number as SeqRecord objects. The keys are the contig names 
		and sequences can be accessed by .seq

		Attributes created using methods:

		GC = percentage of GC content in the genome bin
		bin_length_bp = the total length of the DNA sequence in the
		bin (across all contigs mapped to it)
		number_contigs = the number of contigs mapped to each bin
		N50 = a statistic that can be thought of as the point of
		half of the mass distribution of sequence. Re. the shortest 
		sequence length which covers 50% of the sequence in the genome bin

		*Later Expansion:
		redundancy
		completness
		best_bin
		"""

		self.bin_number = bin_number
		self.contigs = contigs 
		self.GC = None
		self.bin_length_bp = None
		self.number_contigs = None
		self.N50 = None
		self.redundancy = None
		self.completness = None
		self.best_bin = True

	def calculate_bin_descriptors(self):

		"""
		This method calculates the values of the following attributes:

		self.GC
		self.bin_length_bp
		self.number_contigs
		self.N50
		"""

		GC_list = []
		for contig in self.contigs:
			GC_list.append(GC(self.contigs[contig].seq))
		self.GC = sum(GC_list)/len(GC_list)

		self.number_contigs = len(self.contigs)

		contig_lengths = []
		for contig in self.contigs:
			contig_lengths.append(len(self.contigs[contig].seq))
		self.bin_length_bp = sum(contig_lengths)

		sorted_contig_lengths = sorted(contig_lengths)[::-1]
		half_length = (self.bin_length_bp/2)
		lengths = []
		current_length = 0
		N50 = None
		for length in sorted_contig_lengths:
			if sum(lengths) < half_length:
				lengths.append(sorted_contig_lengths[current_length])
				current_length += 1
			else:
				N50 = lengths[-1]
		self.N50 = N50

	def checkM(self):

		"""
		*Later Expansion:
		Will run the program checkM to calculate the redundancy and 
		completness of each bin
		"""

		pass

	def best_bin(self):

		"""
		*Later Expansion:
		Based on user-input (thresholds for completness and redundancy,
		-all flag) will determine is a bin is a best_bin or not 
		(whether to change the True default to False)
		"""

		return self.best_bin

#############################################################################################

#Command-line usage with argparser

parser = argparse.ArgumentParser(prog='find_best_genomebins.py', 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='Generate descriptive information and fasta files for genome bins from \
	an assembly fasta file and associated binning output (contigs mapped to bin numbers)')

parser.add_argument('-assem', '--assembly_path', required=True, help = 'Path to the \
	metagenome assembly file in fasta format with the extension .fna or .fasta')

parser.add_argument('-map', '--contig_to_bin_map_path', required=True, help = 'Path to \
	the mapping file containing contigs with associated bin assignments in the format \
	"contig_name, bin_#", with the "bin_#" as an integer. The contig and bin can be \
	either comma or tab seperated. It should be a text file with either a .txt extension \
	or no extension')

parser.add_argument('-o', '--output_folder', required=True, help = 'Path to the output \
	folder, will be created if it does not already exist')

#parser.add_argument('--compl', type=float, default=0.7, help = 'Completness threshold \
#	a genome bin has to pass to be output and described. As a proportion (float from and \
#	including 0 to 1)')

#parser.add_argument('--redun', type=float, default=1.2, help = 'Redundancy threshold \
#	a genome bin has to pass to be output and described. As a proportion (float >= 1)')

#parser.add_argument('--all', action='store_true', help = 'Include this flag if you want \
#	all bins output and described instead of only those that pass the given completness and \
#	redundancy thresholds (the default)')

#parser.add_argument('--t', type=int, default=1, help = 'Number of threads for checkM \
#	to use when calculated the completness and redundancy of each bin')

args = parser.parse_args()

print "The path given for the assembly: " + args.assembly_path
print "The path given for the mapping file: " + args.contig_to_bin_map_path
#print "The completness threshold is " + args.compl
#print "The redundancy threshold is: " + args.redun
#print "All bins will be output: " + args.all
#print "Number of threads used by CheckM: " + args.t
print "Files will be output in this directory: " + args.output_folder
#parser.print_help()

#############################################################################################

#Acutal implementation

#with open('log.txt', 'w') as log:
#	log.write("Testing \n")
#with open('log.txt', 'a') as log:
#	log.write("1,2,3 \n")

#Initiate the Assembly class
sequences = Assembly(args.assembly_path, args.contig_to_bin_map_path, args.output_folder)

#Parse the assembly
sequences.parse_fasta_file()

#Parse the binning file
sequences.parse_binning_file()

#Initiate the Bin class with bins parsed from the assembly and binning file
sequences.make_Bin_instances()

#Calculate characteristics of each Bin instance
for genome_bin in sequences.bins:
	genome_bin.calculate_bin_descriptors()

#Output table of bin characteristics
sequences.table_bin_characteristics()

#Output plot of bin characteristics
sequences.plot_bin_characteristics()

#Extract genome bins (folder for each bin containing a fasta file of contigs)
sequences.extract_bins()

#log.close()












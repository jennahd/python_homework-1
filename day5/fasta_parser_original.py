#WORK IN PROGRESS
#Script almost BUT not done plot and not fixed for style yet...(will finish tonight)
#Script will be adapted to fit PEP8 standards and docstrings need to be added
#!/usr/bin/env python

"""
This script is designed to parse a given fasta file
To test the functionality of this script execute exercise_day5.py
This second script will import the fasta parser
and test its functionality and behaviour
"""

import os
import re
import sys
import matplotlib

class FastaParser(object):
    """fasta_file_path = the path to the fasta file re. object
    This file takes a given fasta file and will parse the content"""
    def __init__(self, fasta_file):
        if not os.path.exists(fasta_file):
            raise IOError("This file does not exist")
        if fasta_file == "":
            raise TypeError("The file is missing!")
        fasta_file = open(fasta_file, 'r')
        self.fasta_file = fasta_file.read().split('>')
        del self.fasta_file[0]
        sequence_dict = {}
        for line in self.fasta_file:
            header = line.split('\n',1)[0]
            sequence = re.sub('\n', '', (line.split('\n',1)[1]))
            sequence_dict[header] = sequence
        self.sequence_dict = sequence_dict
        self.count = len(self.fasta_file)
    def __len__(self):
        return len(self.fasta_file)
    def __getitem__(self, i):
        if type(i) == int:
            if i <= self.count:
                sequence = self.fasta_file[i]
                return re.sub('\n', '', (sequence.split('\n',1)[1]))
            else:
                raise IndexError("You called a sequence index that is above the range of sequences")
        elif type(i) == str:
            if i in self.sequence_dict:
                return self.sequence_dict[i]
            else:
                raise KeyError("You called a sequence header that does not exist")
    def extract_length(self, length):
        sequences = self.sequence_dict
        filtered_seqs = []
        for i in sequences.values():
            if len(i) < length:
                filtered_seqs.append(i)
        return filtered_seqs
    #def length_dist(self, path):
    #    with open(path, 'w') as plot_file:
    #        lengths = []
    #        number_of_each_length =[]
    #        for i in sequence.values():
    #            lengths.append(count(i))





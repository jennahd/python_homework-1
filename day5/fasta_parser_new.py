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
import matplotlib.pyplot as plt


class FastaParser(object):


    """
    fasta_file_path = the path to the fasta file re. object
    This file takes a given fasta file and will parse the content
    """


    def __init__(self, fasta_file):

        """
        The attributes of the FastaParser class object include:
        fasta_file = the file containing fasta sequences, and a list
        of the sequences that can be accessed by index
        sequence_dict = a dictionary of sequences whereby sequences can
        be accessed by using the appropriate key (re. the sequence 
        header) count = the number of sequences in the fasta file
        """

        if not os.path.exists(fasta_file):
            raise IOError("This file does not exist")
        if fasta_file == "":
            raise TypeError("The file is missing!")
        fasta_file = open(fasta_file, 'r')
        self.fasta_file = fasta_file.read().split('>')
        del self.fasta_file[0]

        sequence_dict = {}
        for line in self.fasta_file:
            header = line.split('\n', 1)[0]
            sequence = re.sub('\n', '', (line.split('\n', 1)[1]))
            sequence_dict[header] = sequence
        self.sequence_dict = sequence_dict

        self.count = len(self.fasta_file)
    
    def __len__(self):

        """
        This re-definied python object will allow you to print the 
        length of the fasta file, re. the number of sequences
        """

        return len(self.fasta_file)
    
    def __getitem__(self, i):

        """
        This re-definied python object will allow you to access 
        sequences in your fasta file by calling an index or the header
        of the sequence.
        """

        if type(i) == int:
            if i <= self.count:
                sequence = self.fasta_file[i]
                return re.sub('\n', '', (sequence.split('\n', 1)[1]))
            else:
                raise IndexError("You called a sequence index that is \
                                  above the range of sequences")
        elif type(i) == str:
            if i in self.sequence_dict:
                return self.sequence_dict[i]
            else:
                raise KeyError("You called a sequence header that does \
                                not exist")
    
    def extract_length(self, length):

        """
        This method allows the user to specify a specific sequence 
        length for which they want to extract corresponding sequences
        """

        sequences = self.sequence_dict
        filtered_seqs = []
        for i in sequences.values():
            if len(i) < length:
                filtered_seqs.append(i)
        return filtered_seqs

    def length_dist(self, mypath):

        """
        This method will return the length distribution of the 
        sequences in the input file as a distribution in a pdf file.
        """
        
        sequences = self.sequence_dict
        sequence_lengths = []
        for i in sequences.values():
            sequence_lengths.append(len(i))
        if not os.path.exists(mypath):
            os.makedirs(mypath)
        fig = plt.hist(sequence_lengths, bins=range(0, max(sequence_lengths)))
        plot(mypath)
        fig.savefig(mypath)



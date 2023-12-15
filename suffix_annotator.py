######################
#Import Dependencies
######################

import os
import subprocess
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Phylo
import matplotlib.pyplot as plt
import pandas as pd
import csv
import re

######################
#Define Core Functions
######################

##################################################################################################################################

def load_suffix_mapping(root_dir):
    suffix_mapping = {}
    try:
        with open(root_dir+"/occurrences.csv", 'r') as csv_file:
            reader = csv.DictReader(csv_file)
            for row in reader:
                occurrence_id = row['occurrence_id']
                country = row['country']
                suffix_mapping[occurrence_id] = country

        print("Loaded suffix mapping:", suffix_mapping)  # Add this line to print the loaded mapping
        return suffix_mapping
    except Exception as e:
        print(f"An error occurred while loading suffix mapping: {e}")
        return {}

##################################################################################################################################

def map_suffix(root_dir, suffix_mapping):
    try:
        with open(root_dir + "/seq_master.fasta", 'r') as input_fasta, open(root_dir + "/suffix_added_seqs.fasta", 'w') as output_fasta:
            for line in input_fasta:
                if line.startswith('>'):
                    header = line.strip()[1:]  # Remove the ">" character
                    header_parts = header.split("_")
                    sequence_name = header_parts[0]  # Extract the first part as the sequence name
                    suffix = suffix_mapping.get(sequence_name, '')  # Get the suffix from the mapping
                    revised_header_parts = [sequence_name]

                    if len(header_parts) > 1:
                        revised_header_parts.extend(header_parts[1:])  # Include the remaining parts of the original header

                    if suffix:
                        revised_header_parts.append(suffix)

                    revised_header = ">" + "_".join(revised_header_parts) + "\n"
                    output_fasta.write(revised_header)
                else:
                    output_fasta.write(line)  # Write sequence lines as is

        print("Mapping completed successfully!")
    except Exception as e:
        print(f"An error occurred: {e}")

##################################################################################################################################

#####################
#Define Main Function
#####################

def main():

    print("Welcome to the SUFFIX ADDER")

    proj_dir = input("Enter the path to your project folder. It should contain 'occurrences.csv' and 'seq_master.fasta':  ")

    suffix_map = load_suffix_mapping(proj_dir)

    print(suffix_map)

    map_suffix(proj_dir, suffix_map)

    print("Suffixes Added!")

##################
#Run Main Function
##################

if __name__ == "__main__":
    main()

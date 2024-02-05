import os
import csv
import requests
import pandas as pd
import wget
import gzip
import shutil
from Bio import SeqIO
import matplotlib.pyplot as plt
from requests.exceptions import HTTPError
import sys
import tkinter as tk
from tkinter import filedialog, messagebox

class GINSAClass:

###############################################################################################################################
    def __init__(self, proj_dir=None, species_name=None):
        """
        Initialize the GINSAClass instance.

        Parameters:
        - proj_dir (str): Path to the project directory.
        - species_name (str): Name of the species for analysis.
        """
        if proj_dir is not None and species_name is not None:
            self.proj_dir = proj_dir
            self.species_name = species_name
            self.gen_sp = species_name.split()

            if len(self.gen_sp) == 2:
                self.genus = self.gen_sp[0]
                self.species = self.gen_sp[1]
            elif len(self.gen_sp) == 1:
                self.genus = self.gen_sp[0]
                self.species = self.gen_sp[0]
            else:
                self.genus = None
                self.species = None
        else:
            self.proj_dir = None
            self.species_name = None
            self.genus = None
            self.species = None

###############################################################################################################################
    def search_species_occurrences(self, species_name, limit=300):
        """
        Search for species occurrences using the GBIF API.

        Parameters:
        - species_name (str): Name of the species for which occurrences are searched.
        - limit (int): Maximum number of occurrences to retrieve per API call.

        Returns:
        - occurrence_ids (list): List of occurrence IDs.
        - occurrence_data (dict): Dictionary containing occurrence data.
        """
        base_species_url = "https://api.gbif.org/v1/species/match"
        base_occurrence_url = "https://api.gbif.org/v1/occurrence/search"

        # Step 1: Use species match to get taxonKey
        species_params = {
            "name": species_name,
            "strict": False, 
        }

        try:
            species_response = requests.get(base_species_url, params=species_params)
            species_data = species_response.json()

            if "usageKey" in species_data:
                taxon_key = species_data["usageKey"]
            else:
                print(f"Could not find taxonKey for species: {species_name}")
                return [], None

        except Exception as e:
            print(f"An error occurred while fetching taxonKey: {e}")
            return [], None

        # Step 2: Use taxonKey to get occurrences
        base_occurrence_url = "https://api.gbif.org/v1/occurrence/search"
        occurrence_params = {
            "taxonKey": taxon_key,
            "limit": limit,
            "fields": "key,country,decimalLatitude,decimalLongitude,identifier,eventDate",
        }

        occurrence_ids = []
        occurrence_data = {"results": []}
        offset = 0

        try:
            while True:
                occurrence_params["offset"] = offset
                occurrence_response = requests.get(base_occurrence_url, params=occurrence_params)
                occurrence_json = occurrence_response.json()

                if "results" in occurrence_json:
                    occurrence_data["results"].extend(occurrence_json["results"])
                    occurrence_ids.extend([occurrence["key"] for occurrence in occurrence_json["results"]])

                    if len(occurrence_json["results"]) < limit:
                        break  # Reached the end of results
                    else:
                        offset += limit
                else:
                    break

        except Exception as e:
            print(f"An error occurred while fetching occurrences: {e}")

        return occurrence_ids, occurrence_data


###############################################################################################################################
    def create_folders(self, occurrence_ids, proj_dir):
        """
        Create folders for each occurrence ID in the project directory.

        Parameters:
        - occurrence_ids (list): List of occurrence IDs.
        - proj_dir (str): Path to the project directory.
        """
        if not os.path.exists(proj_dir):
            os.makedirs(proj_dir)

        for occurrence_id in occurrence_ids:
            folder_path = os.path.join(proj_dir, str(occurrence_id))
            os.makedirs(folder_path, exist_ok=True)

###############################################################################################################################
    def generate_csv(self, occurrence_ids, occurrences, proj_dir):
        """
        Generate a master CSV file containing occurrence information.

        Parameters:
        - occurrence_ids (list): List of occurrence IDs.
        - occurrences (dict): Dictionary containing occurrence data.
        - proj_dir (str): Path to the project directory.
        """
        if not os.path.exists(proj_dir):
            os.makedirs(proj_dir)

        with open(os.path.join(proj_dir, 'occurrences.csv'), mode='w', newline='', encoding='utf-8') as csv_file:
            fieldnames = ['occurrence_id', 'country', 'latitude', 'longitude', 'ENA_index', 'eventDate']
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

            writer.writeheader()

            if isinstance(occurrences, dict) and "results" in occurrences:
                occurrences_list = occurrences["results"]
            else:
                print("Unexpected format of occurrences data. Unable to generate CSV.")
                return
            # Find the occurrence with the matching key
            for occurrence_id in occurrence_ids:
                occurrence_info = next((occurrence for occurrence in occurrences_list if isinstance(occurrence, dict) and occurrence.get("key") == occurrence_id), None)

                if occurrence_info:
                    try:
                        occurrence_id = occurrence_info.get("key")
                        latitude = occurrence_info.get("decimalLatitude")
                        longitude = occurrence_info.get("decimalLongitude")
                        ENA_index = occurrence_info.get("identifier")
                        event_date = occurrence_info.get("eventDate", "")

                        country = occurrence_info.get("country", "").replace(',', '') if occurrence_info.get("country") else None

                        if ENA_index is not None:
                            underscore_index = ENA_index.find('_')
                            ENA_index = ENA_index[:underscore_index] if underscore_index != -1 else ENA_index
                        
                        writer.writerow({
                            'occurrence_id': occurrence_id,
                            'country': country,
                            'latitude': latitude,
                            'longitude': longitude,
                            'ENA_index': ENA_index,
                            'eventDate': event_date  # Included upon reviewer request
                        })

                    except ValueError:
                        print(f"Failed to parse JSON for occurrence {occurrence_id}.")
                else:
                    print(f"Occurrence data not found for ID {occurrence_id}.")

        print("CSV file created successfully.")

###############################################################################################################################
    def process_directory(self, proj_dir):
        """
        Extract and delete compressed files in the project directory.

        Parameters:
        - proj_dir (str): Path to the project directory.
        """
        for root, dirs, files in os.walk(proj_dir):
            for file in files:
                if file.endswith(".gz"):
                    gz_file_path = os.path.join(root, file)
                    extract_path = os.path.splitext(gz_file_path)[0]

                    with gzip.open(gz_file_path, 'rb') as f_in, open(extract_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                    os.remove(gz_file_path)
                    print(f"Extracted and deleted: {gz_file_path}")

###############################################################################################################################
    def is_valid_url(self, url):
        """
        Checks if URL is valid.

        Parameters:
        - url (str): URL to be checked.

        Returns:
        - bool: True if the URL is valid, otherwise returns False.
        """
        try:
            response = requests.head(url)
            response.raise_for_status()  # Raise an exception for HTTP errors
            return True
        except HTTPError as http_err:
            if response.status_code == 404:
                return False

        except Exception as err:
            print(f"An error occurred while checking the URL {url}: {err}")
            return False

###############################################################################################################################
    def ssu_fasta_grab(self, csv_file, proj_dir):
        """
        Download SSU FASTA files from ENA via GBIF for each occurrence found.

        Parameters:
        - csv_file (str): Path to the CSV file containing occurrence metadata.
        - proj_dir (str): Path to the project directory.
        """
        # Load the master CSV file into a Pandas DataFrame
        df = pd.read_csv(csv_file)

        # Replace with ENA API endpoint
        api_base_url = "https://www.ebi.ac.uk/metagenomics/api/v1"

        # Iterate through the DataFrame and process each occurrence
        for ENA_index, occurrence_id in zip(df['ENA_index'], df['occurrence_id']):
            if pd.isna(ENA_index):  # Check if ENA_index is NaN
                print(f"ENA link not found in occurrence metadata for ID {occurrence_id}.")
                continue 

            url = f"{api_base_url}/analyses/{ENA_index}/downloads"

            # Check if the URL is valid
            if not self.is_valid_url(url):
                print(f"Skipping processing for invalid URL: {url}")
                continue

            try:
                response = requests.get(url)
                response.raise_for_status()
                data = response.json()

                print(f"Processing {ENA_index}...")

                # Create directory with name as occurrence_id in project directory
                subdir_path = os.path.join(proj_dir, str(occurrence_id))
                if not os.path.exists(subdir_path):
                    os.makedirs(subdir_path)

                fasta_links = []

                for entry in data['data']:
                    link_entry = entry['links']['self']
                    if link_entry.endswith("SSU.fasta.gz"):
                        fasta_links.append(link_entry)

                if fasta_links:
                    print(f"Downloading {len(fasta_links)} file(s) for {ENA_index}...")
                    for fasta_link in fasta_links:
                        file_name = os.path.basename(fasta_link)
                        file_path = os.path.join(subdir_path, file_name)
                        wget.download(fasta_link, file_path)
                    print(" Download complete.")
                else:
                    print(f"No fasta files found for {ENA_index}. Saving 'no_fasta.txt'...")
                    no_fasta_file = os.path.join(subdir_path, "no_fasta.txt")
                    with open(no_fasta_file, "w") as f:
                        f.write("No fasta files found.")

            except HTTPError as http_err:
                if response.status_code == 404:
                    print(f"HTTP Error 404: ENA data not found for ID {occurrence_id}")
                else:
                    print(f"HTTP Error {response.status_code} occurred: {http_err}")
            except Exception as err:
                print(f"An error occurred while processing ID {occurrence_id}: {err}")
                
###############################################################################################################################
    def mapseq_grab(self, csv_file, proj_dir):
        """
        Download MAPSeq files from ENA via GBIF for each occurrence.

        Parameters:
        - csv_file (str): Path to the CSV file containing occurrence metadata.
        - proj_dir (str): Path to the project directory.
        """
        # Load the master CSV file into a Pandas DataFrame
        df = pd.read_csv(csv_file)

        # Replace with ENA API endpoint
        api_base_url = "https://www.ebi.ac.uk/metagenomics/api/v1"

        # Iterate through the DataFrame and process each occurrence
        for ENA_index, occurrence_id in zip(df['ENA_index'], df['occurrence_id']):
            if pd.isna(ENA_index):  # Check if ENA_index is NaN
                print(f"ENA link not found in occurrence metadata for ID {occurrence_id}.")
                continue 

            if "MGY" not in ENA_index.upper():  # Confirm prefix text points to EMBL ENA
                print(f"ENA link not found in occurrence metadata for ID {occurrence_id}.")
                continue 

            url = f"{api_base_url}/analyses/{ENA_index}/downloads"
            
            try:
                response = requests.get(url)
                response.raise_for_status()
                data = response.json()

                print(f"Processing {ENA_index}...")

                # Create a directory named the occurrence_id in project directory if it does not already exist
                subdir_path = os.path.join(proj_dir, str(occurrence_id))
                if not os.path.exists(subdir_path):
                    os.makedirs(subdir_path)

                mapseq_links = []

                # Check if 'data' key is present in the response dictionary
                if 'data' not in data:
                    print("Error: 'data' key not present in the API response.")
                    continue  # Continue to the next iteration

                # Loop through the data to find links ending with "SSU_MAPSeq.mseq.gz"
                for entry in data['data']:
                    link_entry = entry['links']['self']
                    if link_entry.endswith("SSU_MAPSeq.mseq.gz"):
                        mapseq_links.append(link_entry)

                if mapseq_links:
                    print(f"Downloading {len(mapseq_links)} file(s) for {occurrence_id}...")
                    for mapseq_link in mapseq_links:
                        file_name = os.path.basename(mapseq_link)
                        file_path = os.path.join(subdir_path, file_name)
                        wget.download(mapseq_link, file_path)
                    print(" Download complete.")
                else:
                    print(f"No MAPSeq files found for {occurrence_id}. Saving 'no_mapseq.txt'...")
                    no_mapseq_file = os.path.join(subdir_path, "no_mapseq.txt")
                    with open(no_mapseq_file, "w") as f:
                        f.write("No MAPSeq files found.")

            except HTTPError as http_err:
                if response.status_code == 404:
                    print(f"HTTP Error 404: ENA data not found for ID {occurrence_id}")
                else:
                    print(f"HTTP Error {response.status_code} occurred: {http_err}")
            except Exception as err:
                print(f"An error occurred while processing ID {occurrence_id}: {err}")

###############################################################################################################################
    def find_target_in_mapseq(self, subdirectory, file_path, genus, species):
        """
        Search for target strings in MAPSeq files and save truncated strings.

        Parameters:
        - subdirectory (str): Path to the subdirectory containing MAPSeq files.
        - file_path (str): Path to the MAPSeq file.
        - genus (str): Genus of target species.
        - species (str): Target species.

        Returns:
        - truncated_strings (list): List of truncated strings.
        """
        truncated_strings = [] 

        try:
            with open(file_path, 'r') as file:
                content = file.readlines()

                for line in content:
                    line = line.strip()
                    if genus in line and species in line:
                        # Find the index of the search word in the line
                        word_index = line.index(species)

                        # Extract the preceding string
                        preceding_string = line[:word_index].strip()
                        # Truncate to only the text before the first space
                        truncated_string = preceding_string.split()[0]
                        truncated_strings.append(truncated_string)

            with open(subdirectory+"/truncated_strings.txt", 'w') as output_file:
                for item in truncated_strings:
                    output_file.write("%s\n" % item)

        except FileNotFoundError:
            print("MAPSeq Search Error: The file does not exist or the path is incorrect.")

        return truncated_strings

###############################################################################################################################
    def sift_fasta(self, fasta_file_path, truncated_strings_path):
        """
        Extract DNA sequences from FASTA files based on truncated strings.

        Parameters:
        - fasta_file_path (str): Path to the FASTA file.
        - truncated_strings_path (str): Path to the file containing truncated strings.

        Returns:
        - sifted_data (dict): Dictionary containing truncated strings and corresponding DNA sequences.
        """
        sifted_data = {}  # Dictionary to store truncated_strings and their corresponding DNA sequences

        if truncated_strings_path:
            try:
                with open(truncated_strings_path, 'r') as file:
                    lines = file.readlines()
                lines = [line.strip() for line in lines]
            except FileNotFoundError:
                print("Truncated strings file not found. Continuing without filtering.")
                lines = []

        # Parse the FASTA file and extract sequences
        with open(fasta_file_path, "rt") as fasta_file:
            records = SeqIO.parse(fasta_file, "fasta")
            for record in records:
                sequence_name = record.id
                sequence = str(record.seq)

        # Search for sequence names in the truncated_strings list
                for truncated_string in lines:
                    if truncated_string == sequence_name:
                        sifted_data.setdefault(truncated_string, []).append(sequence)
                        break

        return sifted_data


#######################################################################################################
    def save_sifted_data_to_csv(self, sifted_data, output_file):
        """
        Save sifted data to CSV file.

        Parameters:
        - sifted_data (dict): Dictionary containing truncated strings and DNA sequences.
        - output_file (str): Path to the output CSV file.
        """
        with open(output_file, "w", newline="") as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(["Truncated String", "Sequences"])

            for truncated_string, sequences in sifted_data.items():
                csv_writer.writerow([truncated_string, "\n".join(sequences)])

###############################################################################################################################
    def analyze_subdir_mapseq(self, proj_dir, genus, species):
        """
        Analyze MAPSeq files in subdirectories for target strings.

        Parameters:
        - proj_dir (str): Path to the project directory.
        - genus (str): Genus of target species.
        - species (str): Target species.
        """
        # Iterate through subdirectories in proj_dir
        for subdirectory in os.listdir(proj_dir):
            subdir_path = os.path.join(proj_dir, subdirectory)

            # Skip if not a directory
            if not os.path.isdir(subdir_path):
                print(f"Skipping non-directory: {subdirectory}")
                continue

            print(f"Analyzing subdirectory {subdirectory} for MAPSeq files")

            # Check for required files and text file
            for file_name in os.listdir(subdir_path):
                if file_name.endswith("SSU_MAPSeq.mseq"):
                    print(file_name)

                    self.find_target_in_mapseq(subdir_path, subdir_path+"/"+file_name, genus, species)
                    print("Truncated strings obtained from MAPSeq file.")

###############################################################################################################################
    def analyze_subdir_fasta(self, proj_dir):
        """
        Analyze FASTA files in subdirectories for target strings.

        Parameters:
        - proj_dir (str): Path to the project directory.
        """
        # Iterate through subdirectories in proj_dir
        for subdirectory in os.listdir(proj_dir):
            subdir_path = os.path.join(proj_dir, subdirectory)

            print("subdirectory path:", subdir_path)

            # Skip if not a directory
            if not os.path.isdir(subdir_path):
                print(f"Skipping non-directory: {subdirectory}")
                continue

            print(f"Analyzing subdirectory {subdirectory} for FASTA files.")

            fasta_files = [file_name for file_name in os.listdir(subdir_path) if file_name.endswith("SSU.fasta")]

            # Check if any FASTA files are found in the subdirectory
            if fasta_files:
                # Iterate through FASTA files in the subdirectory
                for fasta_file in fasta_files:
                    fasta_filepath = os.path.join(subdir_path, fasta_file)
                    print("Found FASTA file:", fasta_filepath)

                    # Make path for truncated string file
                    truncated_strings_filepath = os.path.join(subdir_path, "truncated_strings.txt")

                    # Perform analysis using sift_fasta and save_sifted_data_to_csv
                    sifted_data = self.sift_fasta(fasta_filepath, truncated_strings_filepath)

                    # Create a new CSV filename based on the original FASTA filename
                    csv_filename = os.path.splitext(fasta_file)[0] + "_extracted_sequences.csv"

                    # Save the sifted data to the new CSV file
                    self.save_sifted_data_to_csv(sifted_data, os.path.join(subdir_path, csv_filename))

                print("Extraction and CSV creation complete for all FASTA files.")
            else:
                print("No FASTA files found in the subdirectory.")

###############################################################################################################################
    def check_for_csv(self, subdirectory_path):
        """
        Check if CSV files exist in subdirectory.

        Parameters:
        - subdirectory_path (str): Path to the subdirectory.

        Returns:
        - bool: True if CSV files exist, False otherwise.
        """
        csv_files = [file for file in os.listdir(subdirectory_path) if file.endswith('.csv')]
        return len(csv_files) > 0

##################################################################################################################################
    def check_dir(self, proj_dir):
        """
        Check the presence of CSV files in subdirectories and generate a bar plot.

        Parameters:
        - proj_dir (str): Path to the project directory.
        """
        subdirectories = [d for d in os.listdir(proj_dir) if os.path.isdir(os.path.join(proj_dir, d))]

        presence_data = []
        for subdir in subdirectories:
            has_csv = self.check_for_csv(os.path.join(proj_dir, subdir))
            presence_data.append(has_csv)

        num_have_csv = presence_data.count(True)
        num_not_have_csv = presence_data.count(False)

        plt.figure(figsize=(6, 6))
        plt.bar(['With FASTA', 'Without FASTA'], [num_have_csv, num_not_have_csv], color=['green', 'red'])
        #plt.xlabel('Project Sub-Folders') # Un-comment if X-axis label desired
        plt.ylabel('Number of Sub-folders')
        plt.title('Sub-folders Containing Sequences')
        plt.tight_layout()
        plt.savefig(proj_dir + "/sifting_results", dpi=600)

##################################################################################################################################
    def combine_csv_files(self, proj_dir):
        """
        Combine CSV files from subdirectories into a master FASTA file.

        Parameters:
        - proj_dir (str): Path to the project directory.

        Returns:
        - fasta_path (str): Path to the combined FASTA file.
        """
        subdirectories = [d for d in os.listdir(proj_dir) if os.path.isdir(os.path.join(proj_dir, d))]

        combined_sequences = {}

        for subdir in subdirectories:
            if self.check_for_csv(os.path.join(proj_dir, subdir)):
                subdir_sequences = []
                for csv_file in os.listdir(os.path.join(proj_dir, subdir)):
                    if csv_file.endswith('.csv'):
                        csv_path = os.path.join(proj_dir, subdir, csv_file)
                        csv_content = pd.read_csv(csv_path)
                        if 'Sequences' in csv_content.columns:
                            sequences_data = csv_content['Sequences']
                            subdir_sequences.extend(sequences_data)

                if subdir_sequences:
                    combined_sequences[subdir] = subdir_sequences

        if combined_sequences:
            fasta_content = ""
            for subdir, sequences in combined_sequences.items():
                if len(sequences) == 1:
                    fasta_content += f">{subdir}\n{sequences[0]}\n"
                else:
                    for i, sequence in enumerate(sequences, start=1):
                        fasta_content += f">{subdir}_{i}\n{sequence}\n"

            fasta_path = os.path.join(proj_dir, 'seq_master.fasta')
            with open(fasta_path, 'w') as fasta_file:
                fasta_file.write(fasta_content)
            print("Combined sequences saved as 'seq_master.fasta' in the root directory.")

            return fasta_path

        else:
            print("No sequences found in subdirectories.")

##################################################################################################################################
    def seq_master_lengths(self, fasta_file, proj_dir):
        """
        Analyze sequence lengths in a master FASTA file and generate a bar plot.

        Parameters:
        - fasta_file (str): Path to the master FASTA file.
        - proj_dir (str): Path to the project directory.
        """
        # Check if the FASTA file exists
        if not os.path.exists(fasta_file):
            print(f"Error: The file '{fasta_file}' does not exist.")
            return

        try:
            # Read sequences from the FASTA file
            sequences = list(SeqIO.parse(fasta_file, "fasta"))

            # Get sequence names and lengths
            seq_names = [seq.id for seq in sequences]
            seq_lengths = [len(seq) for seq in sequences]

            # Create a bar plot to visualize sequence lengths
            plt.figure(figsize=(15, 9))
            plt.bar(seq_names, seq_lengths, color='blue')
            plt.xlabel('Sequence Names')
            plt.ylabel('Sequence Length')
            plt.title('Sequence Lengths from FASTA File')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            
            # Save the plot to the project directory
            plt.savefig(os.path.join(proj_dir, "sequence_lengths.png"), dpi=600)

        except Exception as e:
            print(f"Error reading sequences from '{fasta_file}': {e}")

##################################################################################################################################
    def analyze_nucleotide_freqs(self, fasta_file, proj_dir):
        """
        Analyze nucleotide frequencies in FASTA file.

        Parameters:
        - fasta_file (str): Path to the FASTA file.
        - proj_dir (str): Path to the project directory.
        """
        sequences = list(SeqIO.parse(fasta_file, "fasta"))

        total_sequence = "".join(str(seq_record.seq) for seq_record in sequences)

        # Calculate nucleotide frequencies for the entire fasta file
        nucleotide_freq = {
            'A': total_sequence.count('A'),
            'T': total_sequence.count('T'),
            'C': total_sequence.count('C'),
            'G': total_sequence.count('G')
        }

        # Define colors for each nucleotide
        colors = ['blue', 'orange', 'green', 'red']

        # Plot nucleotide frequencies with colored bars
        plt.figure(figsize=(8, 6))
        bars = plt.bar(nucleotide_freq.keys(), nucleotide_freq.values(), color=colors)
        plt.xlabel('Nucleotides')
        plt.ylabel('Frequency')
        plt.title('Overall Nucleotide Frequencies')
        plt.tight_layout()

        total_count = sum(nucleotide_freq.values())
        print("Overall Nucleotide Frequencies:")
        for nucleotide, count in nucleotide_freq.items():
            print(f"{nucleotide}: {count} ({(count / total_count) * 100:.2f}%)")
        plt.savefig(proj_dir + "/nucleotide_frequencies", dpi=600)

##################################################################################################################################
    def main(self, proj_dir=None, species_name=None): 
        # Arguments proj_dir and species_name set to None for GUI class interface.
        # Command line version does not include these arguments
        """
        Main function that drives the GINSA analysis process.

        Parameters:
        - proj_dir (str): Path to the project directory.
        - species_name (str): Name of species for analysis.
        """
        if proj_dir is not None and species_name is not None:
            # Use proj_dir and species_name in analysis logic
            self.proj_dir = proj_dir
            self.species_name = species_name
            self.gen_sp = species_name.split()

            if len(self.gen_sp) == 2:
                self.genus = self.gen_sp[0]
                self.species = self.gen_sp[1]
            elif len(self.gen_sp) == 1:
                self.genus = self.gen_sp[0]
                self.species = self.gen_sp[0]
            else:
                self.genus = None
                self.species = None
        else:
            if self.proj_dir is None or self.species_name is None:
                print("Error: project directory and species name are required.")
                return

        gen_sp = species_name.split()

        if len(gen_sp) == 2:
            genus = gen_sp[0]
            print("Genus: ", genus)
            species = gen_sp[1]
            print("Species: ", species)

        elif len(gen_sp) == 1:
            genus = gen_sp[0]
            print("Genus: ", genus)
            species = gen_sp[0]
            print("Only Genus provided. Search will be based on: ", species)

        elif len(gen_sp) > 2:
            print("Too many search terms.")

        else:
            print("No search taxon entered.")

        occurrence_ids, occurrence_data = self.search_species_occurrences( species_name)

        num_occurrences = len(occurrence_ids)

        if num_occurrences > 0:
            print(f"Number of occurrences found: {num_occurrences}")

        else:
            print("No occurrences found for this taxon.")

        if num_occurrences > 300:
            print("There seem to be many ocurrences of this taxon in GBIF. Please ensure that you have sufficient storage.")

        self.create_folders(occurrence_ids, proj_dir)
        print("Subdirectory folders created")
        print("Gathering biogeography data on all the occurrences found in GBIF...")

        self.generate_csv(occurrence_ids, occurrence_data, proj_dir)

        print("Occurrence metadata CSV created")
        csv_file = proj_dir+"/occurrences.csv"

        print("Grabbing FASTA and MAPSeq files from ENA via GBIF...")
        self.ssu_fasta_grab(csv_file, proj_dir)
        self.mapseq_grab(csv_file, proj_dir)

        print("Unpacking compressed files...")
        self.process_directory(proj_dir)

        print("Searching MAPSeq files...")
        self.analyze_subdir_mapseq(proj_dir, genus, species)

        print("Searching FASTA files...")
        self.analyze_subdir_fasta(proj_dir)

        print("Sequence acquisition complete! Analyzing your sequences and saving them to a FASTA file.")

        print("Checking sub-folders for new sequences.")
        self.check_dir(proj_dir)

        print("Aggregating sequences into a master csv file.")
        fasta_path = self.combine_csv_files(proj_dir)

        if fasta_path is None:
            print("Exiting program. No sequences found.")
            sys.exit()

        print("Generating a master FASTA file containing your new sequences.")
        self.seq_master_lengths(fasta_path, proj_dir)

        print("Analyzing the sequences in your master FASTA file.")
        self.analyze_nucleotide_freqs(fasta_path, proj_dir)

class GINSAGUI:
        
    """
    GINSAGUI - GUI for GINSA Analysis

    This class defines a simple GUI for performing analysis using GINSA.

    Parameters:
    - master (tk.Tk): The root window for the GUI.

    Note: Requires the GINSAClass.
    """
    def __init__(self, master):
        """
        Initializes the GINSAGUI instance.
        """
        self.master = master
        self.master.title("GINSA GUI")

        self.project_dir_label = tk.Label(master, text="Project Directory:")
        self.project_dir_label.grid(row=0, column=0, padx=10, pady=5, sticky="e")

        self.project_dir_entry = tk.Entry(master, width=50)
        self.project_dir_entry.grid(row=0, column=1, padx=10, pady=5, columnspan=2, sticky="w")

        self.browse_project_dir_button = tk.Button(master, text="Browse", command=self.browse_project_dir)
        self.browse_project_dir_button.grid(row=0, column=3, padx=5, pady=5, sticky="w")

        self.species_name_label = tk.Label(master, text="Species Name:")
        self.species_name_label.grid(row=1, column=0, padx=10, pady=5, sticky="e")

        self.species_name_entry = tk.Entry(master, width=50)
        self.species_name_entry.grid(row=1, column=1, padx=10, pady=5, columnspan=2, sticky="w")

        self.execute_button = tk.Button(master, text="Execute Analysis", command=self.execute_analysis)
        self.execute_button.grid(row=2, column=1, columnspan=2, pady=10)

    def browse_project_dir(self):
        """
        Opens a file dialog to browse and select the project directory.
        Also updates project directory entry in the GUI.
        """
        directory = filedialog.askdirectory()
        if directory:
            self.project_dir_entry.delete(0, tk.END)
            self.project_dir_entry.insert(0, directory)

    def execute_analysis(self):
        """
        Executes the GINSA analysis using the provided project directory and species name.
        Also displays a message box with the analysis result.
        """
        project_dir = self.project_dir_entry.get()
        species_name = self.species_name_entry.get()

        if not project_dir or not species_name:
            messagebox.showerror("Error", "Please enter both project directory and species name.")
            return

        ginsa_instance = GINSAClass()  # Initialize GINSAClass without parameters; default arguments = None
        try:
            ginsa_instance.main(project_dir, species_name)
            messagebox.showinfo("Analysis Complete", "GINSA analysis completed successfully!")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during analysis: {str(e)}")

def main():
    root = tk.Tk()
    app = GINSAGUI(root)
    root.mainloop()
    ginsa_instance = GINSAClass(proj_dir, species_name)
    ginsa_instance.main()

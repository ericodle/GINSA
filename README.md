<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/github_username/repo_name">
    <img src="https://github.com/ericodle/GINSA/blob/main/docs/logo3.jpg" alt="Logo" width="600" height="300">
  
  </a>

<h3 align="center">GINSA: <ins>G</ins>b<ins>I</ins>f <ins>N</ins>ext-gen <ins>S</ins>equence <ins>A</ins>ccumulator</h3>

  <p align="center">
 
  
</div>


<!-- ABOUT THE PROJECT -->
# About GINSA

GINSA is a Python program that exploits GBIF's biogeographic data and sequence links.
Users input their species of interest, then the program gathers data on each GBIF occurrence.
Data includes the GBIF occurrence ID, latitude/longitude, and country of origin -- all saved to a CSV spreadsheet.
GINSA downloads SSU rDNA sequences from each occurrence, combining them into a single FASTA file.

We tested GINSA on multiple taxa across the Tree of Life.

GINSA is a unique tool utilizing the GBIF interface. In the spirit of GBIF, the code for GINSA is openly available, and we encourage community collaboration. Feel free to fork.

<p align="right">(<a href="#top">back to top</a>)</p>

# Prerequisites
## Install Python3
GINSA was developed using the Python programming language. 
Installing Python, specifically newer versions (>=3.10) is therefore recommended.
Using earlier versions of Python may cause dependency errors.
Python3 can be downloaded at: https://www.python.org/downloads/

## For Mac users
Newer versions of Python lack certain SSL (secure sockets layer) certificates on MacOS.
To correct this, perform the following steps in your working environment:
### Step 1: Install certifi
  ```sh
  python3 -m pip3 install --upgrade pip3 
  pip3 install --upgrade certifi
  ```
### Step 2: Configure system SSL certificates
Find the path to your Python distribution's Certificates.command, then enter the following command.
Our system used Python 3.12, therefore the command looked like this:

  ```sh
  sudo /Applications/Python\ 3.12/Install\ Certificates.command
  ```
Your system should now have the required SSL certificates to run GINSA.

# Quick Start

## Install GINSA using PIP
To install GINSA using PIP, simply open a terminal emulator (command prompt, terminal, etc.) and run:

  ```sh
  pip3 install GINSA
  ```
While not necessary, we recommend using GINSA in a virtual environment.

## Option 1: Run GINSA by graphical user interface (GUI)

This is the point-and-click method of interfacing with software.
In your working environment, simply enter:
  ```sh
  GINSA_gui
  ```
From here, simply select your project directory, type in your taxon of interest, and click "Execute Analysis".

## OPTION 2: Run GINSA by command line interface

This interface method allows more experienced users to run GINSA with text commands.
In your working environment, simply enter GINSA_cli followed by the path to your project directory and your taxon of interest in quotation marks.
It should look something like this:

  ```sh
  GINSA_cli /eo/Desktop/test "Lecudina tuzetae"
  ```

# Description of features

When GINSA completes its run, your project folder will contain unique sub-folders for each occurrence as well as several output files. <br />

Below is a description of each output file using the example search taxon term "Lecudina longissima".

- [ ] sifting_results.png

  This image plots a simple histogram showing how many occurrence sub-folders actually contain FASTA and MAPseq files. Some occurrences in GBIF draw from data types other than next-gen sequences, such as human observation or Sanger sequencing. In the future, we may expand the scope of GINSA to grab Sanger SSU sequences from other online databases as such as the European Nucleotide Archive (ENA), the National Center for Biotechnology Information (NCBI), the DNA Data Bank of Japan (DDBJ), and SILVA (Latin for 'forest').

  <img src="https://github.com/ericodle/GINSA/blob/main/docs/sifting_results.png" alt="seq_master.fasta" width="350" height="350">

- [ ] sequence_lengths.png

  This image visualizes the length of each sequence in the seq_master.fasta file. While the average length is around 200 bp, variability in length can be quite broad.

  <img src="https://github.com/ericodle/GINSA/blob/main/docs/sequence_lengths.png" alt="seq_master.fasta" width="450" height="270">


- [ ] nucleotide_frequencies.png

  This image plots the proportion of each nucleotide, A, T, C, and G, as a frequency histogram. Nucleotide proportions serve as an additional quality check on the extracted sequences.

  <img src="https://github.com/ericodle/GINSA/blob/main/docs/nucleotide_frequencies.png" alt="seq_master.fasta" width="320" height="240">
  

- [ ] <img src="https://github.com/ericodle/GINSA/blob/main/docs/occurrences.csv" alt="occurrences.csv" width="350" height="350">

  This file is a basic spreadsheet that can be opened in Excel or any other CSV reader. Data column 1 contains the GBIF occurrence IDs found during the taxon search, while data column 2 contains the country code from which the sample was taken. Data columns 3 and 4 contain latitude and longitude coordinates for each occurrence. Finally, data column 5, labeled "prefix_text", contains the EMBL ENA project title sifted from the GBIF API for each occurrence. "prefix_text" can be replaced into the ENA API link, leading GINSA to the download page for each SSU FASTA and MAPseq file.


- [ ] <img src="https://github.com/ericodle/GINSA/blob/main/docs/seq_master.fasta" alt="seq_master.fasta" width="350" height="350">

  This file is the primary output of GINSA -- a FASTA file containing all the gathered SSU sequences for all occurrences of a species on GBIF. From this FASTA file, researchers can perform downstream phylogenetic analysis and modify the sequence labels with occurrence metadata from occurrences.csv to explore biogeographic patterns. The supplemental script <img src="https://github.com/ericodle/GINSA/blob/main/SUFFIX_ADDER.py (https://github.com/ericodle/GINSA/blob/main/misc/suffix_adder.py)" alt="SUFFIX_ADDER.py" width="350" height="350"> is included as a tool to add the country code to the end of each sequence name. Moreover, users are encouraged to play with the code and modify it to their particular research needs.



# Key Functions

- [ ] search_species_occurrences()

  This function makes the initial GBIF API call, searching the occurrences database for all cases of the user's target species. search_species_occurrences generates a list of all matching occurrence IDs which will be used as a reference by downstream functions.

- [ ] ssu_fasta_grab()

  This function makes API calls to EMBL's metagenomic data repository, ENA, and replaced the stored variable "ENA_index" for each occurrence into the EMBL URL. Then, the function searches JSON tags on the terminal site to find and download all FASTA files containing SSU contigs. 

- [ ] mapseq_grab()

  This function performs the same actions as ssu_fasta_grab but targets associated MAPseq files rather than FASTA files.

- [ ] find_target_in_mapseq()

  This function opens the MAPseq file for each occurrence sub-directory and searches for sequence labels corresponding to the target species. For example, if you searched "Lecudina longissima",  find_target_in_mapseq will search the MAPseq file for all sequence labels containing the strings "Lecudina" and "longissima. If only a single taxon string is provided, such as for the genus-level search "Labyrinthula", then all sequence labels containing "Labyrinthula" will be stored in a list for downstream reference.

- [ ] sift_fasta()

  This function is complementary to find_target_in_mapseq, and utilizes the sequence label list generated previously to extract each SSU sequence contained in the (often very large) FASTA files grabbed from EMBL. The SSU contig FASTA files provided do not contain the names of the taxa they identify, hence the need for the MAPseq file.

- [ ] combine_csv_files()

  Finally, once the SSU sequences are isolated for each distinct occurrence, they are gathered into a single file named "seq_master.fasta". Researchers may then use this extracted metagenomic data in their research project.

<p align="right">(<a href="#top">back to top</a>)</p>

## Citing our work

Our paper provides an overview of the methodology and results used in this repository.

If you find GINSA useful, we ask that you cite the following: coming soon



<!-- LICENSE -->
## License

This project is open-source and is released under the [MIT License](LICENSE). Feel free to use and build upon our work while giving appropriate credit.



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo_name.svg?style=for-the-badge
[contributors-url]: https://github.com/github_username/repo_name/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/github_username/repo_name.svg?style=for-the-badge
[forks-url]: https://github.com/github_username/repo_name/network/members
[stars-shield]: https://img.shields.io/github/stars/github_username/repo_name.svg?style=for-the-badge
[stars-url]: https://github.com/github_username/repo_name/stargazers
[issues-shield]: https://img.shields.io/github/issues/github_username/repo_name.svg?style=for-the-badge
[issues-url]: https://github.com/github_username/repo_name/issues
[license-shield]: https://img.shields.io/github/license/github_username/repo_name.svg?style=for-the-badge
[license-url]: https://github.com/github_username/repo_name/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/linkedin_username
[product-screenshot]: images/screenshot.png

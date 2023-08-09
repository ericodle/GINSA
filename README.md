<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/github_username/repo_name">
    <img src="https://github.com/" alt="Logo" width="400" height="350">
  </a>

<h3 align="center">GINSA: <ins>G</ins>b<ins>I</ins>f <ins>N</ins>ext-gen <ins>S</ins>equence <ins>A</ins>nalyzer</h3>

  <p align="center">
    GINSA is a tool to help biodiversity researchers extract sequence and location data via GBIF.<br /> <br />
    
  As marine protistologists, the impetus for this project arose from obstacles faced during our own research.
  Thus, we created a Python program that exploits GBIF's wealth of biogeography data and sequence links.
  Users simply type their species of interest, then GINSA gathers metadata on each GBIF occurrence.
  Gathered data includes the GBIF occurrence ID, latitude/longitude, and country of origin -- all saved to a simple spreadsheet.
  GINSA then downloads SSU rDNA sequences for each occurrence, combining them into a single FASTA file.<br /><br />

  GINSA has been tested on multiple protist taxa at both the species and genus level.
  Those taxa include Lecudina longissima (26 occurrences), Lecudina tuzetae (309 occurrences), Labyrinthula (2,603 occurrences), and Symbiodinium (10,426 occurrences).
  During the development and testing phases, we discovered issues with the available data that deserve community attentions.
  Specifically, we noticed two challenges: 1) A lack of SSU amplification site uniformity, and 2) unbalanced biogeographic representation.<br /><br />
  
  The first challenge refers to the relative size of the 18S SSU gene (approx. 1,800 bp) compared to the next-gen SSU sequence lengths (100-300 bp) recovered using GINSA.
  Simply, sequences from different sites/teams/countries amplified different regions of the SSU gene. This is a problem when attempting to construct a phylogenetic tree, since the sequences are not comparable.
  The second challenge refers to sampling bias present in the datasets. Taking Lecudina tuzetae as an example, the majority of sequences were obtained from Germany. Only a handful of sequences were from other countries such as Australia, Great Brittain, the United States, and Taiwan. 
  We hope the revelation of these issues lead to community-driven solutions and, ultimately, a more robust and reliable GBIF.<br /><br />

  GINSA is the first tool of its kind to our knowledge, and is anticipated to save researchers countless hours in manual data collection.
  In the spirit of GBIF, the code for GINSA is openly available, and we encourage collaboration in making it better. There are links below to report a bug or request a new feature.
  Moreover, we invite everyone interested to follow along with the instructions provided below.
  
  <a href="https://github.com/github_username/repo_name/issues">Report Bug</a>
  
  <a href="https://github.com/github_username/repo_name/issues">Request Feature</a>
  </p>
</div>


<!-- ABOUT THE PROJECT -->
## GINSA Workflow

First,

Second,

Third

Fourth....


<p align="right">(<a href="#top">back to top</a>)</p>


## Getting Started

In the spirit of simplicity, GINSA can be run with a single file named RUN_GINSA.py.
This file can be downloaded directly, or you can clone the repository by copy-pasting these commands into your command prompt (terminal).

  ```sh
  # Replace "your_folderpath_here" with the actual folder where you want the project to go.
  cd /your_folderpath_here
  git clone git@github.com:ericodle/GINSA.git
  ```

For your convenience, a Jupyter Notebook is also provided via the Google Colab service.
This option is nice if you want to see the code as the program runs. It also allows you to run the program using Google's computers rather than your own (the software runs in an internet browser).
To use the Jupyter Notebook version of GINSA, simply click the link:** and then click the "Open in Colab" button. *insert image*
You can watch our instructional video for a guided walkthrough of how to use the Jupyter Notebook.

### Install dependencies using pip
Whether you downloaded RUN_GINSA.py or cloned the repository, you will need Python installed on your computer. GINSA was tested on Python3.
Once Python is installed, you will need to copy-paste these commands into your command prompt (terminal).
GINSA uses 

  ```sh
  # Install required packages for GINSA using pip.
  pip install requests
  pip install wget
  pip install biopython
  ```

### Execute RUN_GINSA.py

Create an empty folder (directory) and take note of its full folder path (should look something like "/home/researcher/Desktop/GINSA_run"). All the project files will be saved there.
Ensure you have at least 50 GB of available storage. The actual about of required space will depend on the number of occurrences available for the species you search.
Ensure your internet connection is stable. You will be downloading a lot of large FASTA files. For instance, Lecudina tuzetae had 309 occurrences in GBIF, and RUN_GINSA took about 1 hour to complete.
Ensure you know the full filepath of your RUN_GINSA.py file. We find it easiest to just leave it on the desktop.
Open up a command line prompt (terminal) and execute the following command:

```sh
# Run the Python script titled RUN_GINSA.py
python3 path_to_your_file/RUN_GINSA.py
```

If everything went well, you should be greeted with the welcome message:

```sh
Welcome to the GbIf Next-gen Sequence Analyzer! (GINSA)
```

Followed by a prompt asking for your project directory path:

```sh
Enter the path to an empty folder for this project:/home/researcher/Desktop/lecudina_longissima
```

Next, GINSA will ask for the species you wish to search. 
Here, you can search either the full genus+species name, as in "Lecudina longissima", or just a genus, as in "Labyrinthula".

```sh
Enter your target genus and species (i.e. Lecudina longissima):Lecudina longissima
```

GINSA will respond with some confirmation details:

```sh
Genus:  Lecudina
Species:  longissima
Number of occurrences found: 26
```

GINSA will then proceed with the analysis, informing you of each step being performed in real time.
Note: Analysis time depends on internet connection speed and the number of occurrences in GBIF.

## Completion of Analysis

When GINSA completes its process, you will find that your project folder now contains a unique sub-folder for each occurrence found, as well as several output files. <br />

Below is a description of each output file, using the search term "Lecudina longissima" as an example.

- [ ] occurrences.csv
<img src="https://github.com/" alt="Logo" width="400" height="350">

  What is it?


- [ ] seq_master.fasta
<img src="https://github.com/" alt="Logo" width="400" height="350">

  What is it?


- [ ] nucleotide_frequencies.png
<img src="https://github.com/" alt="Logo" width="400" height="350">

  What is it?


- [ ] sequence_lengths.png
<img src="https://github.com/" alt="Logo" width="400" height="350">

  What is it?


- [ ] sifting_results.png
<img src="https://github.com/" alt="Logo" width="400" height="350">

  What is it?




## Functions

- [ ] train_model.py

This script can be called to train a pre-defined neural network class on labeled MFCC data. Upon training completion, the user will be provided a graph of both training and validation following each train epoch. This graph can be useful in diagnosing neural network issues such as overfitting.

- [ ] test_model.py

This script can be called to test a trained neural network on labeled MFCC data. Once executed, a confusion matrix image will be generated to show the user how well the neural network was able to classify each musical genre.

- [ ] models.py

This bit of code defines the artifical neural network architectures used in our study. Classes for MLP, CNN, LSTM, BiLSTM, and GRU models are written for PyTorch, which we chose over Keras for its greater granular control. Users are welcome to use these model classes to conduct their own experimentation.

- [ ] MFCC_extraction.py

This script extracts MFCCs from the GTZAN dataset music files and saves them in JSON format. The resulting file is about 640 MB in size, and contains an ordered list of 13 MFCC values per segment of each song within the dataset. Moreover, this data is labeled with values 0 through 9 corresponding to one of the ten genres present.

- [ ] MFCC_primer.ipynb

This Jupyter Notebook with direct link to a ready-to-use Google Colab environment is intended to answer the questions "What is an MFCC?" and "How did we get our MFCCs?"

- [ ] live_runs

This folder contains a collection of Jupyter Notebooks (linked to Google Colab) that were saved and uploaded to GitHub immediately after running. In their currently saved state, they serve as a record of the experimental runs on which we base our results. Of course, users are welcome to play around with these scripts and try to beat our top test accuracy of 90.7%!

- [ ] GRU_CM

We herein provide a representative graphic for this README file and, by extension, project. This CM (confusion matrix) obtained upon testing our GRU model serves the additional function of showing our best achieved performance.

- blues.00000.wav

This .wav file is provided as a test music clip. The file itself was taken from the same download of the GTZAN dataset as used throughout this project, and the excerpt is from the classic John Lee Hooker song "One Bourbon, One Scotch, One Beer." It is classified in GTZAN as Blues.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

Contributions make the open source community great. Everyone has a unique combination of skills and experience. Your input is **highly valued**.
If you have ideas for improvement, please fork the repo and create a pull request. 
If this is your first pull request, just follow the steps below:

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the GNU Lesser General Public License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>


Citing
------


EO: lead developer
SR: testing, presentation
KK: testing, image drafting
KW: theory oversight

Manuscript citation coming soon!


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

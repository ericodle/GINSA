<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/github_username/repo_name">
    <img src="https://github.com/" alt="Logo" width="400" height="350">
  </a>

<h3 align="center">GINSA: <ins>G</ins>b<ins>I</ins>f <ins>N</ins>ext-gen <ins>S</ins>equence <ins>A</ins>nalyzer</h3>

  <p align="center">
  GINSA is a tool to help biodiversity researchers extract beiogeographic and taxonomic data via GBIF and associated databases.
    Users simply type their species of interest, then GINSA gathers metadata on each GBIF occurrence of that species.
    Collected metadata includes the GBIF occurrence ID, latitude/longitude, country of collection, and a link to <ins>occurence-specific next-generation sequence repositories</ins>.
    <br />
    <br />
    <a href="https://github.com/github_username/repo_name/issues">Report Bug</a>
    ·
    <a href="https://github.com/github_username/repo_name/issues">Request Feature</a>
  </p>
</div>


<!-- ABOUT THE PROJECT -->
## About this Project

This was my first project with Professor Rebecca Lin (Feng Chia University, Taiwan) during the Taiwan Experience Exchange Program. Classifying music genre was my initial experience both in writing/training ANNs from scratch as well as in audio signal analysis, and I had a lot of catching up to do. Through a wealth of online resources, particularly the MFCC tutorials provided by [Valerio Velardo](https://github.com/musikalkemist), we were able to train models with decent generaliation and test classification accuracy. 

Our results are based on the [GTZAN](http://marsyas.info/index.html) music genre dataset, which provides 10 human-classified genre folders: blues, classical, country, disco, hip-hop, jazz, metal, pop, reggae, and rock. Each genre folder contains 100 30-second audio clips of genre-specific songs in .wav format. Following previous work in this field, we extracted [Mel-frequency cepstral coefficients](https://en.wikipedia.org/wiki/Mel-frequency_cepstrum), or MFCCs, from each audio clip, divided the entire shuffled set into an 80:10:10 train/validation/test split, and played around with multi-layer perceptron, convolutional, and recurrent networks/hyperparamteters until we got a model that achieved at least 90% accuracy.

We hope this project inspires you to contribute to our project, incorporate our tools, and play around with ANN models yourself! 


<p align="right">(<a href="#top">back to top</a>)</p>


## Getting Started

Download this repository by going up to the green "Code" button at the top right and clicking "Download ZIP".

Alternatively, you can also clone the repo directly using the following command.

  ```sh
  # Replace "your_folderpath_here" with the actual folder where you want the project to go.
  cd /your_folderpath_here
  git clone git@github.com:ericodle/GRU_Classifying_GTZAN.git
  ```

> __For this example, the working directory is the repository root directory.__ 

### Install dependencies using pip

  ```sh
  # Install dependencies if necessary. 
  # You may want to work in a virtual environment. Conda environments are nice for that.
  pip install librosa
  pip install torch torchvision
  ```

### Download GTZAN and extract MFCCs

> The training/testing music used in this project comes from the GTZAN music genre dataset, which can be downloaded [here](https://www.kaggle.com/datasets/andradaolteanu/gtzan-dataset-music-genre-classification/download) from Kaggle. 
> You will need to enter some account login details per Kaggle's requirements before the download can begin. 
> Then, you must manually relocate the full dataset (parent folder containing 10 genre sub-folders, each with 100 music clips) into the "GTZAN_dataset" project folder.

```sh
# This script will extract MFCC's from each song clip.
./MFCC_extraction.py
```
> Note that the resulting JSON file is saved in the "MFCCs" folder as a JSON file about 640 MB in size.


### Train a model from scratch

Run the following script to set up a new experiment with default settings.
You will need to specify the type of neural network you want to use.
After the training process is complete, a train/validation curve will be saved in the project root directory.
The final model state will also be saved for the next phase: testing. 

   ```sh
   # Set up a new training run
   ./train_model.py
   ```
Note #1: Training requires a GPU to complete in a timely manner. You can either use your own hardware, or work on a Colab environment.
If you use a GPU, make sure you have cuda and all related dependencies set up in your environment.

Note #2: Training is as much an art as it is a science, and often involves playing around with different hyperparameters. Users are encouraged to go into the train_model.py script and change the optimizer, learning rate, epochs, or other parameters. The default settings represent what worked best for us at the time of experimentation.

### Testing a trained model

You now have a model trained from scratch on MFCCs extracted from the GTZAN music genre dataset. Nice! It is time to see how well it can classify musical genre.
In our conference paper, we used a shuffled 80:10:10 split for training, train phase validation, and testing. Therefore, the music clip segments reserved for testing come from same dataset, but have never been seen by the trained model before. Given the scope of the GTZAN dataset, your trained model is unlikely to distinguish Bunun polyphonic chant music from Ainu rimse dance music. A neural network is only as good as the data on which it is trained. Within the GTZAN training data, how well can your model classify musical genre?

  ```sh
  # Test a pre-trained model.
  ./test_model.py
  ```

Note: The entire MFCC extract JSON file is re-shuffled and split into 80:10:10 train/validation/test subsets each time the train_model.py and test_model.py  scripts are run. Therefore, each train and test run may yield slightly different results. In our experience working on this project, the only factors signifcantly affecting performance were neural network architecture and training hyperparameters.

## Repository Files

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

Please cite the following paper if you use the code provided in this repository.

Conference citation coming soon!


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

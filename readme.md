# STRUCTURE Implementation - CS4775 Final Project

# Overview

This repo contains a partial Python implementation of the [STRUCTURE](http://pritchardlab.stanford.edu/structure.html) algorithm, completed as a required final project for CS4775, taken at Cornell University in Fall 2016. The primary motivation for this project was to implement a non-trivial approach to a problem in the field of computational biology. The objectives were to implement an approach discussed in the literature and apply it to a real dataset.


The implementation can be found in `structureJonathan.py`. The repository also contains the source code for the C implementation of STRUCTURE, which can be downloaded [here](http://pritchardlab.stanford.edu/structure_software/release_versions/v2.3.4/html/structure.html). A PDF writeup of my results can be found in `structureWriteup.pdf`. 

# STRUCTURE Overview
In the field of population genetics, classification of individuals into populations is a common method used to study various problems. STRUCTURE takes a model-based clustering method for using multilocus genotype data to infer structure and assign individuals to populations. More information can be found in the writeup or on the official STRUCTURE website [here](http://pritchardlab.stanford.edu/structure.html). 

# Data
The algorithms are considered on a portion of the Hapmap3 dataset, as acquired by [Admixture](https://www.genetics.ucla.edu/software/admixture/). 


# Some Implementation Details

At the time of completion, the model without admixture was fully implemented, showing good results with both STRUCTURE and [Admixture](https://www.genetics.ucla.edu/software/admixture/). The model with admixture is mostly implemented, but does not show results that agree with STRUCTURE or Admixture. I estimate that the implementation is 90% complete. 


# Statement of Originality

I, Jonathan, am the sole contributor to both the implementation and data analysis. Guidance was generously given by Professor Williams and Melissa Hubisz. 
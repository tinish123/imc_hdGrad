# Code & Data for Computing High-Degree Polynomial Gradients In-Memory

This repo contains the codes and data required to produce results in [1].

[1] Bhattacharya, T., et al. "Computing High-Degree Polynomial Gradients in Memory." arXiv preprint arXiv:2401.16204 (2024).

## Getting Started

### Installation

    git clone https://github.com/tinish123/imc_hdGrad.git

### Dataset

Code to download the remaining 3-SAT instances:

    cd sat
    wget https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/RND3SAT/uf20-91.tar.gz
    wget https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/RND3SAT/uf50-218.tar.gz
    wget https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/RND3SAT/uf75-325.tar.gz
    wget https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/RND3SAT/uf100-430.tar.gz
    tar xvf *.gz
    cd ..

### Description

- The ***sat*** directory contains all the 3-SAT instances that are required for simulations. The randomly generated 14-variable 64-clause 3-SAT instances used for Fig. 5 and Fig. S21 of Supplementary Information, are already provided in the sat directory. All other required instances for N = 20, 50, 75 and 100, need to be downloaded.
- ***data*** directory where all the break and gain-value error models (generated from SPICE simulations pertaining to Fig. S17-S18 of Supplementary File), TTS and run-length data (generated from error-model incorporated WalkSAT/SKC and HO-HNN simulations pertaining to Fig. 5 of main text and FIg. S21 of Supplementary File) and custom polynomial file for real-valued gradient computation demonstration are stored.
- ***new_data*** sub-directory where all new simulation result files can be directed to be stored in.
- *walksat_package.py*, package for performing WalkSAT/SKC simulations with or without error models of break-values and generating and saving uniform random k-SAT problem instances.
- *walksat_experiments.ipynb*, notebook containing experiments related to WalkSAT/SKC simulations. The experiments scipted here can be used to generate all the TTS, run-length data that is used by the other notebook file to make plots.
- *generate_SAT.ipynb*, notebook to generate single or multiple k-SAT instances of arbitrary configurations.
- *hw_model_package.py*, package containing energy and area data of the different circuit components of WalkSAT/SKC and High-Order HNN (HO-HNN) solver, (extracted from Tables S1 and S3 of Supplementary File).
- *TTS_ETS_Plots.ipynb* notebook to generate the final plots corresponding to Fig. 5 of main text and Fig. S21 of Supplementary Information. It also contains code to use the hardware modeling functions in the hw_model_package.py file to obtain energy, power, area results for the WalkSAT/SKC, High-Order and Second-Order HNN solvers, reported in Table S5.
- *real_grad_package.py* contains all relevant fuctions and classes required for real-valued gradient computation demonstration. The constructor function of the class 'real_grad', contains all parameters listed in Table S6 of Supplementary File.
- *real_valued_gradient.ipynb* that contains (i) sample code to run the model for the real-valued gradient computing hardware for a single variable vector, (ii) code to generate Fig. S24 of Supplementary File and (iii) code to explore tradeoff between monomial coefficient precision and maximum order (K) that can be implemented.
- The *custom4_4.poly* file in the ***data*** directory is a compact representation of the polynomial considered in Fig. S23 of Supplementary File. First row is the header and the last two numbers denote the values of N and M respectively. Each subsequent row corresponds to a monomial. All numbers before the '0' represent the indices of variables that are part of the product of the monomial and the number right after the '0' represents the coefficient. The last zero of the file is to be ignored.


### Contact

Tinish Bhattacharya (tinish@ucsb.edu, mymail.tinishbhat@gmail.com)

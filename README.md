# Code & Data for Computing High-Degree Polynomial Gradients In-Memory

This repo contains the codes and data required to produce results in [1].

[1] Bhattacharya, T., et al. "Computing High-Degree Polynomial Gradients in Memory." arXiv preprint arXiv:2401.16204 (2024).

## Getting Started

### Installation

    git clone https://github.com/tinish123/imc_hdGrad.git

### Dataset

The randomly generated 14-variable 64-clause 3-SAT instances used for Fig. 5 and Fig. S21 of Supplementary Information, are already provided in the sat directory. All other required instances for N = 20, 50, 75 and 100, need to be downloaded:

    cd sat
    wget https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/RND3SAT/uf20-91.tar.gz
    wget https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/RND3SAT/uf50-218.tar.gz
    wget https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/RND3SAT/uf75-325.tar.gz
    wget https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/RND3SAT/uf100-430.tar.gz
    tar xvf *.gz
    cd ..

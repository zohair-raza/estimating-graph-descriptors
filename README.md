# Estimating Descriptors for Large Graphs

This repository includes the code for the paper "Estimating Descriptors for Large Graphs" - accepted for publication at the 24th Pacific-Asia Conference on Knowledge Discovery and Data Mining (arXiv link: https://arxiv.org/abs/2001.10301).

The base code was written by Kijung Shin for their work on Tri-Fly, and retrieved from https://github.com/kijungs/trifly. Please read their help files if you have trouble running the code.

GABE - Estimates a graph descriptor based on the Graphlet Kernel

MAEVE - Estimate a graph descriptor based on NetSimile

## Running the code

To run either program, go to their root directory, run 'make', and run the following command:

mpirun -n [#number of processes] ./bin/mpi  --trials [#trials] --budget [#number of edges to store on each worker] [edge list file] [output directory]

The output is a directory with two files: one file has the time taken, and the other has the graph embedding.

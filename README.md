# ATLAS Anti-Kt Algorithm Testing
This project was built to implement the Anti-kt algorithm for jet reclustering, along with different reclustering metrics, for the ATLAS experiment in C++. The goal was to write its function in a low-level with the intention to later implement it in on an FPGA. The Anti-kt algorithm was made to conform with the implementation in the [fastjet library](https://fastjet.fr)'s [documentation](https://fastjet.fr/repo/fastjet-doc-3.4.3.pdf).

## Requirements and Setup
The project runs using CERN's [ROOT library](https://root.cern.ch/). For any file you wish to run, make sure to have the command line running them be currently in the same directory (to ensure relative file reading has the right folder). For the folder labelled "3_TestJetMathAndAlgos", which tests this algorithm implementation against the fastjet library's, you must have fastjet installed and then compile fastjet_example.cpp using the command at the top of the file. Follow the guides on the following pages to install these libraries:

1. [ROOT install](https://root.cern.ch/install/)
2. [fastjet install](https://fastjet.fr/quickstart.html)

## Description
The main implementation of the Anti-kt algorithm is in JetMath.h in the top directory.

TODO

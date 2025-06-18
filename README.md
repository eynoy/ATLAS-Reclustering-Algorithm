# ATLAS Anti-Kt Algorithm Testing
This project was built to implement the Anti-kt algorithm for jet reclustering, along with different reclustering metrics, for the ATLAS experiment in C++. The goal was to write its function in a low-level with the intention to later implement it on an FPGA (i.e. in VHDL/Verilog). The Anti-kt algorithm was made to conform with the implementation in the [fastjet library](https://fastjet.fr)'s [documentation](https://fastjet.fr/repo/fastjet-doc-3.4.3.pdf).

## Requirements and Setup
The project runs using CERN's [ROOT library](https://root.cern.ch/). For any file you wish to run, make sure to have the command line running them be currently in the same directory (to ensure relative file reading has the right folder). Follow the guides on the following pages to install these libraries:

1. [ROOT install](https://root.cern.ch/install/)
2. [fastjet install](https://fastjet.fr/quickstart.html)

For the folder labelled "3_TestJetMathAndAlgos", which tests this algorithm implementation against the fastjet library's, you must have fastjet installed and then compile fastjet_example.cpp using the command at the top of the file.

## Description
The main implementation of the Anti-kt algorithm is in JetMath.h in the top directory.

1. The first folder, "**1_CreateDataSets**", is the algorithm for creating a "Data Set" custom file from ROOT tree nodes, for replicable testing. The file format consists of:

| Item | Number of Bytes |
| ----------- | ----------- |
|  Number of Jets (```n```) | ```sizeof(int)``` |
| Array of Eta Values | ```n*sizeof(double)``` |
| Array of Phi Values | ```n*sizeof(double)``` |
| Array of Transverse Momentum Values | ```n*sizeof(double)``` |
| Array of Mass Values | ```n*sizeof(double)``` |

2. The second folder, "**2_VisulaizeDataSet**", is the algorithm for visualizing the jets in a "Data Set" file. It includes options for a histogram of the mass values, a 2D plot showing where jets are in eta-phi space, and an incomplete 3D visulaization of the jets.

3. The third folder, "**3_TestJetMathAndAlgos**", compares three versions of the Anti-kt algorithm against each other to measure that they agree and their relative performance. The first is the implementation in the fastjet library, the second is an initial implementation of the algorithm, and the third is a more optimized implementation that stores metric values of pairs of jets between passes.

4. The fourth folder, "**4_CreateRocCurves**", finally tests the algorithm according to a chosen reclustering metric and generates a [ROC](https://en.wikipedia.org/wiki/Receiver_operating_characteristic) comparison of a True Positive Rate (TPR) of keeping desired jet clusters against a False Positive Rate of keeping "background" jet clusters, simulated by appropriate data sets.

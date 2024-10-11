# Efficient Approximation of Centrality Measures in Uncertain Graphs
## Bachelors Thesis by Daniel Ketels, May 2023

This repository contains the implementation of both the PSP-harmonic and PSP-betweenness heuristics, as well as the implementation of Monte Carlo for harmonic and betweenness on uncertain graphs. It was designed by Daniel Ketels for his Bachelors Thesis on the efficient approximation of betweenness and harmonic centrality in uncertain graphs. 

The code base consists of three .cpp source files and two .hpp header files. Additionally, there is a makefile to compile a binary 'UncertainCentrality'. In the subdirectory 'testcases', there is a python script 'test.py' and a set of 30 end to end testcases. The PSP algorithms are tested using hard coded results (as the results are deterministic). Monte Carlo is tested by printing the sampled graphs to the respective output files, parsing them into NetworKit and then using the NetworKit centraliy scores as reference. 

The code generally assumes valid input. While it does handle some invalid parameters and graph files, I did not spend much time on this part and it is surely not complete, neither did I test invalid input.

### main.cpp
The 'main.cpp' is responsible for parsing the input parameters, calling the algorithms and outputting the result. It takes 4 or 5 arguments from the commandline in this order:

    - 1. the hyperparameter phi (float)
    - 2. the amount of Monte Carlo samples (uint64_t)
    - 3. (optional) one of the distance functions: d_er (default), d_med, d_maj
    - 4. the algorithm to run: harmonic or betweenness (or distance, which is a function only used in the testsuite)
    - 5. a *.txt file describing an uncertain graph

The input file *.txt needs to have the following format:

    n
    v1,u1,p1
    v2,u2,p2
    ...
    vk,uk,pk

n (uint64_t) is the amount of nodes in the graph; vi, ui (both uint64_t) are node indices with (with 0<=vi,ui<n). They are connected by an uncertain edge with probability pi (double) with 0<=pi<=1.

Additionally, it contains 5 compile time paramers, which are set by "#define X" at the top of the 'main.cpp':

    - OUTPUT_PRECISION (uint): fixed number of decimal places in cenrality score output
    - OUTPUT_MC_SAMPLES (bool): if set to true, the sampled MC graphs are writtne to the output (just for the testsuite)
    - SAVE_MEMORY (bool): call slower PSP version that deallocate memory if possible
    - STDOUT (bool): if set to true, write result to stdout (for simexpal), else there will be a result file
    - SHELL_OUTPUT: if this is defined, the current status of the algorithm is traced in the shell

Note that the testsuite only works properly if we have OUTPUT_PRECISION = 5, OUTPUT_MC_SAMPLES = true, STDOUT = false. 

After parsing all input paramets, the main function calls the selected algorithms and then prints the formated result. It also measures the execuation time of each algorithm. 

To create the data structure for an uncertain graph, it includes the 'UncertainGraph.hpp' header. To execute said algorithms, it includes the 'UncertainGraphAlgorithms.hpp' header.

### UncertainGraph
The header file 'UncertainGraph.hpp' contains the interface of the data structure for uncertain graphs, which is implemented in the 'UncertainGraph.cpp' file. The functionality of this data structure only contains what is neccesarry for the algorithms that are considered here (and it would need to be extended on to give a general purpose data structure for uncertain graphs).

### UncertainGraphAlgorithms
After creating an object that represents an uncertain graph, the 'main.cpp' calls functions declaed in the 'UncertainGraphAlgorithms.hpp' header. Only the baseline MonteCarlo functions and baseline PSP-heuristic are accessed by the main.cpp file. However, a variety of different functions are additionally declared in the 'UncertainGraphAlgorithms.hpp' header and implemented in the 'UncertainGraphAlgorithms.cpp' file. This was done to achieve more clarity (readability, maintanability), e.g. by defining the AllShortestPaths algorithms as seperate functions.

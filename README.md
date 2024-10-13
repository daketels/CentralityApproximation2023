# Efficient Approximation of Centrality Measures in Uncertain Graphs
## Code For The Bachelors Thesis by Daniel Ketels, May 2023

This repository contains the C++ & OpenMP implementation of the PSP-harmonic and PSP-betweenness heuristics (Possible Shortest Paths)f for harmonic closeness and betweenness centrality estimates on uncertain graphs. 
I completed this codebase in May 2023, after working on the neede the theoretical framework, espacially since this problwm was not explored when I choose it for my Bachelors Thesis.
The Title of the thesis is "Efficient approximation of Centrelity Measures in uncertain graphs" (graphs where the existence of edges is probabilistic).

Hence, the time devoted to on this C++ project was limited, espacially as the 50 pages of theoretical work where my priority.
Because of a multiple relevamt scientific discoveys in the relatively new field of uncertain graphs, the thesis was made publicily aviable (open source, CC lincences) at multipl places, e.g. in the section of combinatorical mathematics at ArXiv.org: 
https://arxiv.org/abs/2409.18062 Still, while gaining further knowledge in high performance C++ and making first contacts with OpenMP parallelization, 
I took work on this code serously. Anything else would be highly negligent as well, as the experiments performed using this code where 
supposed to evaluate the accuracy and runtime of my nove algorithm for Harmonic as well as for a reevaluation of the PSP-Betweeness experiments (which sadly did not proove to be reproducible).

To use ground truth centralites with high accucarcy, Monte Caro with a big sample size was as exmplyed (also implemented here and parallelized with OpenMP).
The code itself consists of three .cpp source files and two .hpp header files. Additionally, there is a makefile to compile the binary 'UncertainCentrality'.
In the subdirectory 'testcases', there is a python script 'test.py' and a set of 30 end to end testcases. The PSP algorithms where tested using hard coded (exact) results (hence, on smaller graphs).
Descriptions of all 30 testcases (E2E Tests) and TikZ Drwaings of them can be found as well.

The code generally assumes valid input parameters. While it does handle some invalid parameters and graph files, I did not spend much time on this part and it this whole project surely is complete in that regard, 
as I designed everything just as an (efficient) exmpirical "proof" of the algorithmic concepts.

In the main.cpp file, the input parameters (and the graph input file) is handdled, the algorithms are run and the results writte with multiple output options.
In the header UncertainGraph.hpp, an interface for barebone data structure for uncertain graphs is given, which is implemented in UncertainGraphs.cpp.
Most work was in implementing the multititure of different algorithms and their variants, efficient and in parallel, furthermore designe to be tested on a system of distributed computers.
All algorithm implementations can be found in UncertainGraphAlgorithms.cpp, which of course also has a header file UncertainGraphAlgorithms.hpp.

Keep in mind that some algorithms where experimental, some just for testing or debugging purposes and, most importantly, this was not meant to be either published or to part of any cooperation.
Overall, UncertainGraphAlgorithms.cpp contains roughly 1,400 lines of code, while the other 4 files make up about 550 additional lines. 
Furthermore, I wrote multiple python scripts for data conversion, statistical analysis and plotting of experimental results.
If, by any slim chance some one tries to read this code, feel free to contact me on GitHub (or e.e. on ResearchhGate or Google Scholar).

Additionally, I wrote thousands of lines of Python code to generate, convert, statistically valuate and draw plots. However, I don't tihnk this could be of any interest.

### main.cpp paramters
The 'main.cpp' is responsible for parsing the input parameters, calling the algorithms and outputting the result. It takes 4 or 5 arguments from the commandline in this order:

    -  the hyperparameter phi (float)
    -  the amount of Monte Carlo samples (uint64_t)
    -  (optional) one of the distance functions: d_er (default), d_med, d_maj
    -  the algorithm to run: harmonic or betweenness (or distance, which is a function only used in the testsuite)
    -  a *.txt file describing an uncertain graph

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

#include <iostream>
#include <fstream>
#include <chrono>
#include "UncertainGraph.hpp"
#include "UncertainGraphAlgorithms.hpp"

// compile time parameters
//      the testsuite requires:
//          OUTPUT_PRECISION = 5
//          OUTPUT_MC_SAMPLES = true
//          !defined STDOUT
//          !defined SHELL_OUTPUT
//      for simexpal experiments:
//          OUTPUT_PRECISION = 30
//          OUTPUT_MC_SAMPLES = false
//          defined STDOUT
//          !defined SHELL_OUTPUT
#define OUTPUT_PRECISION 5             // fixed number of decimal places in the output
#define OUTPUT_MC_SAMPLES true         // if set to true: write MC samples output file (just for testing purposes)
#define SAVE_MEMORY false              // use slower allShortestPaths algorithms that deallocates memory when possible
#define STDOUT false                   // write result to stdout instead of result file (for simexpal)

//#define SHELL_OUTPUT // if defined: trace current status in shell
#if defined SHELL_OUTPUT
    #define SHELL(msg) std::cout << msg;
#else
    #define SHELL(msg)
#endif

int main(int argc, char *argv[]) {
    if (argc > 6 || argc < 5){
        std::cerr << "Exactly 4 or 5 arguments have to be passed:\n"
                  << "1. phi (non negative decimal) for PSP\n"
                  << "2. mcSampleCount (integer) number of samples for Monte Carlo (can be 0 -> no MC)\n"
                  << "3. (optional) distance function (can be d_er = default, d_maj, d_med)\n"
                  << "4. algorithm to run (can be harmonic, betweenness, distance (for tests only))\n"
                  << "5. path to *.txt describing an Uncertain Graph (relative from executable)\n\n"
                  << "format of file:\n"
                  << "---------------------------------\n"
                  << "nodeCount (positive int)\n"
                  << "first edge\n"
                  << ".\n"
                  << ".\n"
                  << ".\n"
                  << "last edge\n"
                  << "---------------------------------\n"
                  << "where edge = \"int,int,float\"  and always\n"
                  << "0<=int<nodeCount (node indices), 0<=float<=1 (edge prob)\n";
        return 1;
    }
    // hasNoDist = 1 -> use d_er (default), no distance function was provided
    const int hasNoDist = argc == 5 ? 1 : 0;

    double phi;
    std::stringstream ss(argv[1]);
    ss >> phi;
    if (ss.fail() ||phi < 0.0){
        std::cerr << "failed to read valid phi\n";
        return 1;
    }

    uint64_t mcSampleCount;
    ss = std::stringstream(argv[2]);
    ss >> mcSampleCount;
    if (ss.fail()){
        std::cerr << "failed to read valid MC sample size\n";
        return 1;
    }

    std::string argString;

    // set function pointer to selected distance function (d_er: default)
    double (*distFunc)(const std::vector<double>&);
    distFunc = getDistanceER;
    if (!hasNoDist){
        argString = argv[3];
        if (argString == "d_maj"){
            distFunc = getDistanceMaj;
        }
        else if (argString == "d_med") {
            distFunc = getDistanceMed;
        }
        else if (argString != "d_er"){
            std::cerr << "invalid distance function was provided: " << argString << "\n";
            return 1;
        }
    }

    enum algorithm {Harmonic, Betweenness, Distance};
    algorithm alg;
    argString = argv[4 - hasNoDist];
    if (argString == "harmonic"){
        alg = Harmonic;
    }
    else if(argString == "betweenness"){
        alg = Betweenness;
    }
    else if (argString == "distance"){
        alg = Distance;
    }
    else{
        std::cerr << "invalid algorithm was provided: " << argString << "\n";
        return 1;
    }

    // read and construct graph (i.e. read input file line by line, then call edge list constructor))
    SHELL("trying to read graph from file: " << argv[5 - hasNoDist] << "...\n")
    std::ifstream inFile;
    inFile.open(argv[5 - hasNoDist]);
    if (!inFile.is_open()){
        std::cerr << "Unable to open graph file\n";
        return 1;
    }

    // read node count
    std::string currentLine;
    if (!std::getline(inFile, currentLine)){
        std::cerr << "invalid graph file\n";
        return 1;
    }
    uint64_t nodeCount;
    ss = std::stringstream(currentLine);
    ss >> nodeCount;
    if (ss.fail()){
        std::cerr << "invalid graph file\n";
        return 1;
    }
    if (nodeCount < 2){
        std::cerr << "graph needs to have at least two nodes (3 for betweenness)\n";
        return 1;
    }
    else if (nodeCount < 3 && alg == Betweenness){
        std::cerr << "graph needs to have at least three nodes for betweenness\n";
        return 1;
    }

    //remaining (non-empty) lines contain one edge each of kind 'uint,uint,float'
    edge currentEdge{};
    std::vector<edge> edges{};
    while (std::getline(inFile, currentLine) && !currentLine.empty()){

        uint64_t firstComma, secondComma;
        if ( (firstComma = currentLine.find(',')) == std::string::npos){
            std::cerr << "invalid graph file\n";
            return 1;
        }
        if ( (secondComma = currentLine.find(',', firstComma + 1)) == std::string::npos){
            std::cerr << "invalid graph file\n";
            return 1;
        }

        ss = std::stringstream(currentLine.substr(0, firstComma));
        ss >> currentEdge.v;
        if (ss.fail()){
            std::cerr << "invalid graph file\n";
            return 1;
        }

        ss = std::stringstream(currentLine.substr(firstComma + 1, secondComma - firstComma - 1));
        ss >> currentEdge.u;
        if (ss.fail()){
            std::cerr << "invalid graph file\n";
            return 1;
        }

        ss = std::stringstream(currentLine.substr(secondComma + 1));
        ss >> currentEdge.probability;
        if (ss.fail()){
            std::cerr << "invalid graph file\n";
            return 1;
        }

        if (currentEdge.u == currentEdge.v){
            continue;
        }
        edges.emplace_back(currentEdge);
    }

    inFile.close();
    UncertainGraph G(nodeCount, edges);
    SHELL("done reading and constructing graph with n = " << G.getNodeCount() << ", m = " << G.getEdgeCount() << "...\n")

    std::vector<double> pspResult;
    std::vector<double> mcResult;
    std::chrono::milliseconds durPsp;
    std::chrono::milliseconds durMc;

    // just used for the testsuite, else they remain empty
    std::string mcSampleString;
    std::string outputDistance;

    switch (alg){
        case Harmonic : {
            SHELL("running PSP harmonic...\n")
            auto start = std::chrono::high_resolution_clock::now();
            pspResult = SAVE_MEMORY  ?
                        parallelHarmonicPSPMemory(G, phi, distFunc)
                        : parallelHarmonicPSP(G, phi, distFunc);
            auto end = std::chrono::high_resolution_clock::now();
            durPsp = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            SHELL("PSP done after " << durPsp.count() << " milliseconds...\n")

            if (mcSampleCount == 0){
                SHELL("MC skipped (zero samples provided)\n")
                durMc = std::chrono::milliseconds::zero();
                break;
            }

            SHELL("running MC harmonic...\n")
            start = std::chrono::high_resolution_clock::now();
            mcResult = OUTPUT_MC_SAMPLES ?
                       parallelMcHarmonicPrintSamples(G, mcSampleCount, mcSampleString)
                       : parallelMcHarmonic(G, mcSampleCount);
            end = std::chrono::high_resolution_clock::now();
            durMc = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
            SHELL("MC done after " << durMc.count() << " milliseconds...\n")
        } break;
        case Distance : {
            // just for the testsuite
            outputDistance = distanceTest(G, phi, OUTPUT_PRECISION);
        } break;
        case (Betweenness) : {
            SHELL("running PSP betweenness...\n")
            auto start = std::chrono::high_resolution_clock::now();
            pspResult = SAVE_MEMORY ?
                        parallelBetweennessPSPMemory(G, phi)
                        : parallelBetweennessPSP(G, phi);
            auto end = std::chrono::high_resolution_clock::now();
            durPsp = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            SHELL("PSP done after " << durPsp.count() << " milliseconds...\n")

            if (mcSampleCount == 0){
                SHELL("MC skipped (zero samples provided)\n")
                durMc = std::chrono::milliseconds::zero();
                break;
            }

            SHELL("running MC betweenness...\n")
            start = std::chrono::high_resolution_clock::now();
            mcResult = OUTPUT_MC_SAMPLES ?
                       parallelMcBetweennessPrintSamples(G, mcSampleCount, mcSampleString)
                       : parallelMcBetweenness(G, mcSampleCount);
            end = std::chrono::high_resolution_clock::now();
            durMc = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
            SHELL("MC done after " << durMc.count() << " milliseconds...\n")
        } break;
    }

    // if STDOUT is defined as true: write result to stdout (useful for simexpal)
    // else, for input file *.txt, create file *_result.txt and write output there

    std::ofstream outFile;
    if (alg == Distance || !STDOUT){
        argString = argv[5 - hasNoDist];
        argString = argString.substr(0, argString.find_last_of('.')) + "_result.txt";
        outFile.open(argString);
        if (!outFile.is_open()){
            std::cerr << "unable to open result file: " << argString << "\n";
            return 1;
        }
        if (alg == Distance){
            outFile << outputDistance;
            return 0;
        }
    }
    std::ostream& outStream = STDOUT ? std::cout : outFile;

    ss = std::stringstream();
    ss.precision(OUTPUT_PRECISION);
    ss << std::fixed << pspResult[0];
    outStream << ss.str();
    for (uint64_t i = 1; i < pspResult.size(); ++i){
        ss = std::stringstream();
        ss.precision(OUTPUT_PRECISION);
        ss << std::fixed << pspResult[i];
        outStream << "," << ss.str();
    }
    outStream << "\n";

    if(!mcResult.empty()){
        ss = std::stringstream();
        ss.precision(OUTPUT_PRECISION);
        ss << std::fixed << mcResult[0];
        outStream << ss.str();
        for (uint64_t i = 1; i < mcResult.size(); ++i){
            ss = std::stringstream();
            ss.precision(OUTPUT_PRECISION);
            ss << std::fixed << mcResult[i];
            outStream << "," << ss.str();
        }
    }
    outStream << "\n";
    outStream << "durPSP: " << durPsp.count() << "\n";
    outStream << "durMC: " << durMc.count() << "\n";

    if constexpr (OUTPUT_MC_SAMPLES){
        outStream << mcSampleString;
    }
    return 0;
}
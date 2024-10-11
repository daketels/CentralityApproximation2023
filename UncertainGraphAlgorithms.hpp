#ifndef UNCERTAINCENTRALITY_UNCERTAINGRAPHALGORITHMS_HPP
#define UNCERTAINCENTRALITY_UNCERTAINGRAPHALGORITHMS_HPP

#include "UncertainGraph.hpp"
#include <unordered_set>

/*
from UncertainGraph.hpp:
    using node = uint64_t;
    using nodePair = std::pair<node, node>;
    using nodePairHash = boost::hash<nodePair>;
 */

// explores all shortest source->target paths in G (ignoring probabilities) without using edges in hash set edgeRemoved
// returns tuple of kind:
//      std::vector<double>:    list of absolute existence probabilities of all shortest paths, size = amount of paths
//      uint64_t:               length of the shortest paths
//      std::vector<nodePair>:  vector of edges to remove, one with minimal probability per retrieved path, contains no duplicates
std::tuple<std::vector<double>, uint64_t, std::vector<nodePair>> allSpIgnoringUncertaintyDistance(
        const UncertainGraph &G, node source, node target,
        const std::unordered_set<nodePair, nodePairHash> &edgeRemoved);

// slower version of allSpIgnoringUncertaintyDistance(...) that deallocates memory during the BFS (dequeue depth k+1 -> free mem on every node with depth k)
std::tuple<std::vector<double>, uint64_t, std::vector<nodePair>> allSpIgnoringUncertaintyDistanceMemory(
        const UncertainGraph &G, node source, node target,
        const std::unordered_set<nodePair, nodePairHash> &edgeRemoved);

// calculates vector of G.nodeCount doubles (estimated distance probability distribution p_st)
// hyperparameter phi: PSP stops exploring s,t paths once phi_st (estimated connection probability of s,t) is at least phi
std::vector<double> pspDistanceDistribution (const UncertainGraph &G, node s, node t, double phi);

// version of pspDistanceDistribution(...) with slower (but memory saving) exploration algorithm allSpIgnoringUncertaintyDistanceMemory(...)
std::vector<double> pspDistanceDistributionMemory (const UncertainGraph &G, node s, node t, double phi);

// calculates expected reliable distance between nodes s, t based on given estimated probability distribution p_st
// takes std::vector<double> D where D[i] is the estimate of p_st(i+1), i=0,...,|V|-2 and D[|V|-1] is the estimate of p_st(inf)
double getDistanceER(const std::vector<double> &distribution);

// calculates majority distance between nodes s, t based on given estimated probability distribution p_st
// takes std::vector<double> D where D[i] is the estimate of p_st(i+1), i=0,...,|V|-2 and D[|V|-1] is the estimate of p_st(inf)
double getDistanceMaj(const std::vector<double> &distribution);

// calculates median distance between nodes s, t based on given estimated probability distribution p_st
// takes reference to std::vector<double> D where D[i] is the estimate of p_st(i+1), i=0,...,|V|-2 and D[|V|-1] is the estimate of p_st(inf)
double getDistanceMed(const std::vector<double> &distribution);

// calculates (non normalized) harmonic closeness centrality of node s based on distances in dist
// takes reference to vector<vector<node>> dist where for nodes s < t: dist[s][t-s-1] = d(s,t)
double getHarmonic(node s, uint64_t nodeCount, const std::vector<std::vector<double>> &dist);

// estimates (normalized) harmonic closeness centrality of all nodes in G using distFunc d_er / d_med / d_maj based on PSP in parallel
// hyperparameter phi: PSP stops exploring s,t paths once phi_st (estimated connection probability of s,t) is at least phi
std::vector<double> parallelHarmonicPSP(const UncertainGraph &G, double phi, double (*distFunc)(const std::vector<double>&));

// version of parallelHarmonicPSP(...) that uses slower (but memory saving) exploration algorithm allSpIgnoringUncertaintyDistanceMemory(...)
std::vector<double> parallelHarmonicPSPMemory(const UncertainGraph &G, double phi, double (*distFunc)(const std::vector<double>&));

// traverses G using BFS from source node, calculates vector dist where for each node v: dist[v] = d(source, v)
// only edges in set sampledEdges are used
std::vector<uint64_t> allDistancesSingleSource(const UncertainGraph &G, node source, const std::unordered_set<nodePair, nodePairHash>& sampledEdges);

// calculates harmonic closeness of all nodes in G using parallelized monte carlo
std::vector<double> parallelMcHarmonic(const UncertainGraph& G, uint64_t mcSampleCount);

// calculates harmonic closeness of all nodes in G using parallelized monte carlo
// writes all sampled graphs to mcSampleString (for testing purposes only)
std::vector<double> parallelMcHarmonicPrintSamples(const UncertainGraph& G, uint64_t mcSampleCount, std::string& mcSampleString);

// explores all shortest paths in G (ignoring probabilities) without using edges in set edgeRemoved
// calculates pairs of kind
//     - vector of pairs (Pr(P), B_P) where
//         - Pr(P): absolute probability of path P (double)
//         - B_P: vector<bool> of size G.nodeCount with vec[v] = true if v is an inner node on P
//     - vector of edges to remove, one with minimal probability per retrieved path, contains no duplicates
std::pair<std::vector<std::pair<double, std::vector<bool>>>, std::vector<nodePair>> allSpIgnoringUncertaintyBetweenness(
        const UncertainGraph &G, node source, node target, const std::unordered_set<nodePair, nodePairHash> &edgeRemoved);

// slower version of allSpIgnoringUncertaintyBetweenness(...) that deallocates memory during the BFS (dequeue depth k+1 -> free mem on every node with depth k)
std::pair<std::vector<std::pair<double, std::vector<node>>>, std::vector<nodePair>> allSpIgnoringUncertaintyBetweennessMemory(
        const UncertainGraph &G, node source, node target, const std::unordered_set<nodePair, nodePairHash> &edgeRemoved);

// estimates (normalized) betweenness centrality for all nodes in G using PSP method
// hyperparameter phi: PSP stops exploring more s,t paths once phi_st (estimated connection probability of s,t) is at least phi
std::vector<double> parallelBetweennessPSP(const UncertainGraph &G, double phi);

// version of parallelBetweennessPSP(...) that uses slower (but memory saving) exploration algorithm allSpIgnoringUncertaintyBetweennessMemory(...)
std::vector<double> parallelBetweennessPSPMemory(const UncertainGraph &G, double phi);

// calculates betweenness for all nodes in G using only edges in sampledEdges with brandes algorithm
// the calculated betweenness values are added to the given input vector 'result', reference to it is returned
std::vector<double>& brandesBetweennessAddToAllNodes(const UncertainGraph &G, const std::unordered_set<nodePair, nodePairHash> &sampledEdges, std::vector<double>& result);

// calculates betweenness of all nodes in G using parallelized monte carlo and brandes algorithm
std::vector<double> parallelMcBetweenness(const UncertainGraph& G, uint64_t mcSampleCount);

// calculates betweenness of all nodes in G using parallelized monte carlo and brandes algorithm
// writes all sampled graphs to mcSampleString (for testing purposes only)
std::vector<double> parallelMcBetweennessPrintSamples(const UncertainGraph& G, uint64_t mcSampleCount, std::string& mcSampleString);

// calculates the estimated distance probability function p_st for all s,t in G
// then it uses the p_st values to calculate all three distance functions d_er, d_med, d_maj
// all p_st and distance values are returned in string for testing/verification purposes
std::string distanceTest(UncertainGraph& G, double phi, uint32_t OUTPUT_PRECISION);

#endif //UNCERTAINCENTRALITY_UNCERTAINGRAPHALGORITHMS_HPP
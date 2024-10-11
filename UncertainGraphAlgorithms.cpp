#include "UncertainGraphAlgorithms.hpp"
#include <algorithm>
#include <queue>
#include <numeric>
#include <random>
#include <omp.h>

// Maximal number of threads for OMP
#define OMP_MAX_THREADS omp_get_max_threads()
// "infinity" (maximum double: constexpr)
#define INFINITY_DOUBLE std::numeric_limits<double>::max()
// "infinity" (maximum uint64_t: constexpr)
#define INFINITY_UINT_64 std::numeric_limits<uint64_t>::max()


std::tuple<std::vector<double>, uint64_t, std::vector<nodePair>> allSpIgnoringUncertaintyDistance(
        const UncertainGraph &G, const node source, const node target, const std::unordered_set<nodePair, nodePairHash>& edgeRemoved){

    const uint64_t nodeCount = G.getNodeCount();
    std::vector<node> dist(nodeCount, INFINITY_UINT_64);
    dist[source] = 0;

    // init vector of the shortest path existence probabilities.
    // initially the only shortest path is (source -> source) with probability 1.0; all other nodes hold an empty list of probabilities
    std::vector<std::vector<double>> pathProbabilities(nodeCount);
    pathProbabilities[source].push_back(1.0);

    // minEdges[v] = (edge e, probability p, depth d)
    // at any point, e is an edge with minimal prob. p on all currently known shortest source -> v paths, that was found most recently (at depth d)
    std::vector<std::tuple<nodePair, double, uint64_t>> minEdges(nodeCount, std::make_tuple(nodePair{0, 0}, INFINITY_DOUBLE, 0));

    std::queue<node> Q;
    Q.push(source);
    while (!Q.empty()) {
        const node current = Q.front();
        if (current == target){
            break;
        }
        Q.pop();

        for (const node child : G.neighbors(current)){
            if (edgeRemoved.find(child < current ? nodePair{child, current} : nodePair{current, child}) != edgeRemoved.end()){
                continue;
            }
            if(dist[child] == INFINITY_UINT_64){
                Q.push(child);
                dist[child] = dist[current] + 1;

                // check if {current,child} is minimal prob edge for all (current) shortest source -> child paths
                // else minEdges[child] = minEdges[current]
                const double newProb = G.getProb(current, child);
                const double prevMinProb = std::get<1>(minEdges[current]);
                minEdges[child] = prevMinProb >= newProb ?
                        std::make_tuple(nodePair{current, child}, newProb, dist[child])
                        :  minEdges[current];

                // child inherits all shortest paths from current (i.e. the existence probabilities) but they are multiplied with P(current,child)
                pathProbabilities[child].reserve(pathProbabilities[current].size());
                std::transform(pathProbabilities[current].begin(),
                               pathProbabilities[current].end(),
                               std::back_inserter(pathProbabilities[child]),
                               [newProb](double element) -> double {return element * newProb;});
            }
            else if (dist[child] == dist[current] + 1){
                // check if {current,child} is new minimal edge
                // else check whether the min. edge of current or child should be kept at minEdges[child]
                const double newProb = G.getProb(current, child);
                const double prevMinProbChild = std::get<1>(minEdges[child]);
                const double prevMinProbCurrent = std::get<1>(minEdges[current]);

                if (std::min(prevMinProbChild, prevMinProbCurrent) >= newProb) {
                    minEdges[child] = std::make_tuple(nodePair{current, child}, newProb, dist[child]);
                } else if (prevMinProbChild > prevMinProbCurrent) {
                    minEdges[child] = minEdges[current];
                } else if (prevMinProbChild == prevMinProbCurrent) {
                    // check if dist(source -> min_edge_current) > dist(source -> min_edge_child)
                    if (std::get<2>(minEdges[current]) > std::get<2>(minEdges[child])) {
                        minEdges[child] = minEdges[current];
                    }
                }

                // increase size of pathProbabilities[child] vector to insert all shortest paths that lead over current
                pathProbabilities[child].reserve(pathProbabilities[child].size() + pathProbabilities[current].size());
                // insert their probability multiplied by P(current,child)
                std::transform(
                        pathProbabilities[current].begin(),
                        pathProbabilities[current].end(),
                        std::back_inserter(pathProbabilities[child]),
                        [newProb](double element) -> double {return element * newProb;}
                );
            }
        }
    }
    // source, target disconnected; return empty path prob list, inf dist, empty E_min list
    if (dist[target] == INFINITY_UINT_64){
        return std::make_tuple(std::vector<double>{}, INFINITY_UINT_64, std::vector<nodePair>{});
    }
    // special case where edge {source,target} exists. return ({P(source,target)}, dist = 1, E_min = {{source,target}}); no need for E_min retrieval
    else if (dist[target] == 1) {
        return std::make_tuple(std::move(pathProbabilities[target]), 1, std::vector<nodePair>{nodePair{source, target}});
    }

    // no pre allocation for E_min; generally we will delete much fewer edges than we get paths in big graphs
    std::vector<nodePair> edgesToRemove;
    std::vector<bool> visited(nodeCount, false);
    visited[target] = true;

    // delete {target, child} if this edge is minimal (and not deleted already)
    // else enqueue child, if edge {target,child} is not removed already and child is predecessor on shortest path(s) to source
    Q = std::queue<node>{};
    for (const node child: G.neighbors(target)) {
        if (edgeRemoved.find(target < child ? nodePair{target, child} : nodePair{child, target}) != edgeRemoved.end()){
            continue;
        }
        if (dist[child] == dist[target] - 1) {
            visited[child] = true;
            if (G.getProb(target, child) <= std::get<1>(minEdges[child])) {
                edgesToRemove.emplace_back(target, child);
            } else {
                Q.push(child);
            }
        }
    }
    // traverse shortest paths backwards until min edge is found again; then stop traversal there and add edge to E_min
    while(!Q.empty()){
        const node current = Q.front();
        Q.pop();

        for (const node child : G.neighbors(current)){
            if (edgeRemoved.find(child < current ? nodePair{child, current} : nodePair{current, child}) != edgeRemoved.end()){
                continue;
            }
            if (dist[child] == dist[current] - 1 ){
                const auto [u, v] = std::get<0>(minEdges[current]);
                if (u == child && v == current){
                    edgesToRemove.emplace_back(u, v);
                }
                else if (!visited[child]){
                    visited[child] = true;
                    Q.push(child);
                }
            }
        }
    }
    return std::make_tuple(std::move(pathProbabilities[target]), dist[target], std::move(edgesToRemove));
}

std::tuple<std::vector<double>, uint64_t, std::vector<nodePair>> allSpIgnoringUncertaintyDistanceMemory(
        const UncertainGraph &G, node source, node target, const std::unordered_set<nodePair, nodePairHash> &edgeRemoved){
    // i only comment on the differences to allSpIgnoringUncertaintyDistance(...), more detail is there
    const uint64_t nodeCount = G.getNodeCount();
    std::vector<node> dist(nodeCount, INFINITY_UINT_64);
    dist[source] = 0;

    std::vector<std::vector<double>> pathProbabilities(nodeCount);
    pathProbabilities[source].push_back(1.0);

    // maximal depth of dequeued node yet found, init. zero (depth of source)
    uint64_t currentMaxDepth = 0;
    // vector holds all nodes that have been found at currentMaxDepth (init. just source)
    std::vector<node> nodesWithCurrentDepth;
    nodesWithCurrentDepth.push_back(source);


    std::vector<std::tuple<nodePair, double, uint64_t>> minEdges(nodeCount, std::make_tuple(nodePair{0, 0}, INFINITY_DOUBLE, 0));

    std::queue<node> Q;
    Q.push(source);
    while (!Q.empty()) {
        const node current = Q.front();
        if (current == target){
            break;
        }
        Q.pop();
        // node with depth k+1 is dequeued for the first time; release the path probabilities from all nodes with depth k
        if (dist[current] > currentMaxDepth){
            for (node v : nodesWithCurrentDepth){
                pathProbabilities[v].clear();
                pathProbabilities[v].shrink_to_fit();
            }
            // set new currentMaxDepth and vector nodeWithCurrentDepth
            currentMaxDepth = dist[current];
            nodesWithCurrentDepth.clear();
            nodesWithCurrentDepth.push_back(current);
        }
        for (const node child : G.neighbors(current)){
            if (edgeRemoved.find(child < current ? nodePair{child, current} : nodePair{current, child}) != edgeRemoved.end()){
                continue;
            }
            if(dist[child] == INFINITY_UINT_64){
                Q.push(child);
                dist[child] = dist[current] + 1;

                const double newProb = G.getProb(current, child);
                const double prevMinProb = std::get<1>(minEdges[current]);
                minEdges[child] = prevMinProb >= newProb ?
                        std::make_tuple(nodePair{current, child}, newProb, dist[child])
                        : minEdges[current];

                pathProbabilities[child].reserve(pathProbabilities[current].size());
                std::transform(pathProbabilities[current].begin(),
                               pathProbabilities[current].end(),
                               std::back_inserter(pathProbabilities[child]),
                               [newProb](double element) -> double {return element * newProb;});
            }
            else if (dist[child] == dist[current] + 1){
                const double newProb = G.getProb(current, child);
                const double prevMinProbChild = std::get<1>(minEdges[child]);
                const double prevMinProbCurrent = std::get<1>(minEdges[current]);

                if (std::min(prevMinProbChild, prevMinProbCurrent) >= newProb) {
                    minEdges[child] = std::make_tuple(nodePair{current, child}, newProb, dist[child]);
                } else if (prevMinProbChild > prevMinProbCurrent) {
                    minEdges[child] = minEdges[current];
                } else if (prevMinProbChild == prevMinProbCurrent) {
                    // if dist(source -> min_edge_current) > dist(source -> min_edge_child)
                    if (std::get<2>(minEdges[current]) > std::get<2>(minEdges[child])) {
                        minEdges[child] = minEdges[current];
                    }
                }

                pathProbabilities[child].reserve(pathProbabilities[child].size() + pathProbabilities[current].size());
                std::transform(
                        pathProbabilities[current].begin(),
                        pathProbabilities[current].end(),
                        std::back_inserter(pathProbabilities[child]),
                        [newProb](double element) -> double {return element * newProb;}
                );
            }
        }
    }
    if (dist[target] == INFINITY_UINT_64){
        return std::make_tuple(std::vector<double>{}, INFINITY_UINT_64, std::vector<nodePair>{});
    }
    else if (dist[target] == 1) {
        return std::make_tuple(std::move(pathProbabilities[target]), 1, std::vector<nodePair>{nodePair{source, target}});
    }

    std::vector<nodePair> edgesToRemove;
    std::vector<bool> visited(nodeCount, false);
    visited[target] = true;

    Q = std::queue<node>{};
    for (const node child: G.neighbors(target)) {
        if (edgeRemoved.find(target < child ? nodePair{target, child} : nodePair{child, target}) != edgeRemoved.end()){
            continue;
        }
        if (dist[child] == dist[target] - 1) {
            visited[child] = true;
            if (G.getProb(target, child) <= std::get<1>(minEdges[child])) {
                edgesToRemove.emplace_back(target, child);
            } else {
                Q.push(child);
            }
        }
    }
    while(!Q.empty()){
        const node current = Q.front();
        Q.pop();

        for (const node child : G.neighbors(current)){
            if (edgeRemoved.find(child < current ? nodePair{child, current} : nodePair{current, child}) != edgeRemoved.end()){
                continue;
            }
            if (dist[child] == dist[current] - 1 ){
                const auto [u, v] = std::get<0>(minEdges[current]);
                if (u == child && v == current){
                    edgesToRemove.emplace_back(u, v);
                }
                else if (!visited[child]){
                    visited[child] = true;
                    Q.push(child);
                }
            }
        }
    }
    return std::make_tuple(std::move(pathProbabilities[target]), dist[target], std::move(edgesToRemove));
}

std::vector<double> pspDistanceDistribution (const UncertainGraph &G, const node s, const node t, const double phi){
    const uint64_t nodeCount = G.getNodeCount();

    std::unordered_set<nodePair, nodePairHash> removedEdges;
    removedEdges.reserve(G.getEdgeCount());

    std::vector<double> distribution(nodeCount, 0.0);

    double negatedProbabilityProduct = 1.0;     // prod (1-Pr(s)) over all found the shortest s-t paths s
    double phiST = 0.0;                         // (est.) connection prob. phi_st
    double probabilitySum = 0.0;                // current sum of p_st(k) for k=1,...,|V|-1

    while (phiST < phi){
        // explore next set of paths
        auto [pathProbabilities, dist, edgesToRemove] =
                allSpIgnoringUncertaintyDistance(G, s, t, removedEdges);

        // s,t disconnected
        if (dist == INFINITY_UINT_64){
            break;
        }

        // next p_st probability value: p_st(k) = negatedProbabilityProduct * sum_{S in P_st} Pr(S)
        double newProb = negatedProbabilityProduct *
                         std::accumulate(pathProbabilities.begin(), pathProbabilities.end(), 0.0, std::plus<>());

        // special case where we would get sum p_st(k) >= 1
        if (probabilitySum + newProb >= 1.0){
            distribution[dist - 1] = 1.0 - probabilitySum;
            return distribution;
        }

        // set next p_st value; update p_st sum; update negatedProbabilityProduct
        distribution[dist - 1] = newProb;
        probabilitySum += newProb;
        for (const double prob : pathProbabilities){
            negatedProbabilityProduct *= (1.0 - prob);
        }

        // update phiST; remove E_min edges
        phiST = 1.0 - negatedProbabilityProduct;
        for (const auto &[u, v] : edgesToRemove){
            removedEdges.insert(u < v ? nodePair{u, v} : nodePair{v, u});
        }
    }
    // p_st(inf) = 1 - sum(p_st(k))
    distribution[nodeCount - 1] = 1.0 - probabilitySum;
    return distribution;
}

std::vector<double> pspDistanceDistributionMemory (const UncertainGraph &G, const node s, const node t, const double phi){
    // i only comment on the differences to pspDistanceDistribution(...), more detail is there
    const uint64_t nodeCount = G.getNodeCount();

    // no reserve here; rehash instead if needed, saves memory
    std::unordered_set<nodePair, nodePairHash> removedEdges;
    std::vector<double> distribution(nodeCount, 0.0);

    double negatedProbabilityProduct = 1.0;
    double phiST = 0.0;
    double probabilitySum = 0.0;
    while (phiST < phi){
        // call allSpIgnoringUncertaintyDistanceMemory(...) version of exploration
        auto [pathProbabilities, dist, edgesToRemove] =
                allSpIgnoringUncertaintyDistanceMemory(G, s, t, removedEdges);
        if (dist == INFINITY_UINT_64){
            break;
        }
        double newProb = negatedProbabilityProduct *
                         std::accumulate(pathProbabilities.begin(), pathProbabilities.end(), 0.0, std::plus<>());
        if (probabilitySum + newProb >= 1.0){
            distribution[dist - 1] = 1.0 - probabilitySum;
            return distribution;
        }
        distribution[dist - 1] = newProb;
        probabilitySum += newProb;
        for (const double prob : pathProbabilities){
            negatedProbabilityProduct *= (1.0 - prob);
        }
        phiST = 1.0 - negatedProbabilityProduct;
        for (const auto &[u, v] : edgesToRemove){
            removedEdges.insert(u < v ? nodePair{u, v} : nodePair{v, u});
        }
    }
    distribution[nodeCount - 1] = 1.0 - probabilitySum;
    return distribution;
}

double getDistanceER(const std::vector<double>& distribution){
    // straightforward calculation of exp. rel. distance according to formula
    const uint64_t nodeCount = distribution.size();

    if (distribution[nodeCount - 1] == 1.0){
        return INFINITY_DOUBLE;
    }

    double result = 0.0;
    for (uint64_t i = 0; i < nodeCount - 1; ++i){
        result += distribution[i] * static_cast<double>(i+1);
    }
    result *= (1.0 / (1.0-distribution[nodeCount - 1]));
    return result;
}

double getDistanceMed(const std::vector<double>& distribution){
    // calculation of d_med with small alterations that i thought might be semantically advantageous (but not tested yet)
    const uint64_t nodeCount = distribution.size();

    for (uint64_t dist = 0; dist < nodeCount - 1; ++dist){
        if (distribution[dist] > 0.5){
            return static_cast<double>(dist+1);
        }
    }
    if (distribution[nodeCount - 1] > 0.5){
        return INFINITY_DOUBLE;
    }
    if (distribution[0] == 0.5){
        return 1.0;
    }

    uint64_t dist = 1;
    double sum = distribution[0];
    while (dist < nodeCount){
        sum += distribution[dist];
        if (sum >= 0.5){
            if (sum > 0.5){
                return static_cast<double>(dist);
            }
            break;
        }
        ++dist;
    }
    return static_cast<double>(dist+1);
}

double getDistanceMaj(const std::vector<double>& distribution){
    // simple calculation of d_maj distance according to formula
    const uint64_t nodeCount = distribution.size();

    uint64_t dist = 1;
    double max = distribution[0];
    for (uint64_t i = 1; i < nodeCount; ++i){
        if (max < distribution[i]){
            dist = i + 1;
            max = distribution[i];
        }
    }
    if (dist == nodeCount){
        return INFINITY_DOUBLE;
    }
    return static_cast<double>(dist);
}

double getHarmonic(const node s, const uint64_t nodeCount, const std::vector<std::vector<double>>& dist){
    // dist(s,t) = dist[s][t-s-1] for s < t
    // simply calculates harmonic of node s (using distances in dist) by the usual formula

    double result = 0.0;
    for (node t = 0; t < s; ++t){
        const double distST = dist[t][s - t - 1];
        if (distST != INFINITY_DOUBLE){
            result += 1.0 / distST;
        }
    }
    for (node t = s+1; t < nodeCount; ++t){
        const double distST = dist[s][t - s - 1];
        if (distST != INFINITY_DOUBLE){
            result += 1.0 / distST;
        }
    }
    return result;
}

std::vector<double> parallelHarmonicPSP(const UncertainGraph &G, const double phi, double (*distFunc)(const std::vector<double>&)){
    const uint64_t nodeCount = G.getNodeCount();

    // prepare triangular distance matrix, i.e. |V|-1 vectors of sizes 1, 2, ... , |V| - 1
    std::vector<std::vector<double>> dist;
    dist.resize(nodeCount - 1);
    for (node v = 0; v < nodeCount - 1; ++v){
        dist[v].resize(nodeCount - v - 1);
    }

    // loop over s,t in parallel to get est. p_st distribution and then s-t distance (based on p_st and given distFunc)
#pragma omp parallel for
    for (node s = 0; s < nodeCount - 1; ++s){
        for (node t = s+1; t < nodeCount; ++t){
            std::vector<double> distribution = pspDistanceDistribution(G, s, t, phi);
            dist[s][t-s-1] = distFunc(distribution);
        }
    }

    // loop over all nodes in parallel and calculate their normalized PSP harmonic values based on the distances in dist[][]
    std::vector<double> resultPSP;
    resultPSP.resize(nodeCount);
    const double normalizationFactor = 1.0 / static_cast<double>(nodeCount-1); // divide only once to save time
#pragma omp parallel for
    for (node s = 0; s < nodeCount; ++s){
        resultPSP[s] = normalizationFactor * getHarmonic(s, nodeCount, dist);
    }
    return resultPSP;
}

std::vector<double> parallelHarmonicPSPMemory(const UncertainGraph &G, const double phi, double (*distFunc)(const std::vector<double>&)){
    // the only difference to parallelHarmonicPSP(...) is that we call pspDistanceDistributionMemory(...) when calculating distances
    const uint64_t nodeCount = G.getNodeCount();

    std::vector<std::vector<double>> dist;
    dist.resize(nodeCount - 1);
    for (node v = 0; v < nodeCount - 1; ++v){
        dist[v].resize(nodeCount - v - 1);
    }

#pragma omp parallel for
    for (node s = 0; s < nodeCount - 1; ++s){
        for (node t = s+1; t < nodeCount; ++t){
            std::vector<double> distribution = pspDistanceDistributionMemory(G, s, t, phi);
            dist[s][t-s-1] = distFunc(distribution);
        }
    }

    std::vector<double> resultPSP;
    resultPSP.resize(nodeCount);
    const double normalizationFactor = 1.0 / static_cast<double>(nodeCount-1);
#pragma omp parallel for
    for (node s = 0; s < nodeCount; ++s){
        resultPSP[s] = normalizationFactor * getHarmonic(s, nodeCount, dist);
    }
    return resultPSP;
}

std::vector<uint64_t> allDistancesSingleSource(const UncertainGraph &G, const node source, const std::unordered_set<nodePair, nodePairHash>& sampledEdges){
    // straightforward BFS with check for "e not in sampledEdges ?"
    std::vector<node> dist(G.getNodeCount(), INFINITY_UINT_64);
    dist[source] = 0;

    std::queue<node> Q;
    Q.push(source);

    while (!Q.empty()){
        const node current = Q.front();
        Q.pop();
        for (const node child : G.neighbors(current)){
            if (sampledEdges.find(current < child ? nodePair{current, child} : nodePair{child, current}) == sampledEdges.end()){
                continue;
            }
            if (dist[child] == INFINITY_UINT_64){
                dist[child] = dist[current] + 1;
                Q.push(child);
            }
        }
    }
    return dist;
}

std::vector<double> parallelMcHarmonic(const UncertainGraph& G, const uint64_t mcSampleCount){
    const uint64_t nodeCount = G.getNodeCount();
    const uint64_t edgeCount = G.getEdgeCount();

    // partial harmonic sum results; initially vector of |V| zeros for every thread
    std::vector<std::vector<double>> accumulatedResultsPerThread(OMP_MAX_THREADS, std::vector<double>(nodeCount, 0.0));

    // sample mcSampleCount graphs in parallel, distribute them over all available OMP threads
#pragma omp parallel for num_threads(accumulatedResultsPerThread.size())
    for (uint64_t runNum = 0; runNum < mcSampleCount; ++runNum) {
        const auto threadId = static_cast<std::size_t>(omp_get_thread_num());

        std::random_device rand;
        std::default_random_engine eng(rand());
        std::uniform_real_distribution<double> distribution(0.0, 1.0);

        std::unordered_set<nodePair, nodePairHash> sampledEdges;
        sampledEdges.reserve(edgeCount);

        for (node current = 0; current < nodeCount - 1; ++current) {
            for (const node neighbor: G.neighbors(current)) {
                if (current > neighbor) { // dont sample edges twice
                    continue;
                }
                if (distribution(eng) < G.getProb(current, neighbor)) {
                    sampledEdges.insert(nodePair{current, neighbor});
                }
            }
        }

        // iterate over all nodes v. calculate dist[] with dist[s]=d(s,v) for all s using BFS; calculate H(v) with dist
        for (node v = 0; v < nodeCount; ++v){
            double harmonic = 0.0;
            std::vector<uint64_t> dist = allDistancesSingleSource(G, v, sampledEdges);
            for (node s = 0; s < v; ++s){
                if (dist[s] != INFINITY_UINT_64){
                    harmonic += 1.0 / static_cast<double>(dist[s]);
                }
            }
            for (node s = v+1; s < nodeCount; ++s){
                if (dist[s] != INFINITY_UINT_64){
                    harmonic += 1.0 / static_cast<double>(dist[s]);
                }
            }
            // add resulting harmonic of v to partial sum of current thread
            accumulatedResultsPerThread[threadId][v] += harmonic;
        }
    }
    // iterator starts at second element; first element (partial sum of thread with id 0) is moved into mcResult
    auto it = std::next(accumulatedResultsPerThread.begin());
    std::vector<double> mcResult(std::move(accumulatedResultsPerThread[0]));

    // add partial results of all other threads to mcResult
    while (it != accumulatedResultsPerThread.end()){
        std::transform(it -> begin(), it -> end(),
                       mcResult.begin(), mcResult.begin(), std::plus<>());
        ++it;
    }

    // normalize and get MC average, i.e. divide by (|V|-1)(mcSampleCount)
    const double normalizationAndAveragingFactor = 1.0 / (static_cast<double>((nodeCount - 1) * mcSampleCount));
    std::transform(mcResult.begin(), mcResult.end(), mcResult.begin(),
                   [normalizationAndAveragingFactor](double elem){return elem * normalizationAndAveragingFactor;});
    return mcResult;
}

std::vector<double> parallelMcHarmonicPrintSamples(const UncertainGraph& G, uint64_t mcSampleCount, std::string& mcSampleString){
    // this is just for debugging/testing; same as parallelMcHarmonic(...) with additional output (generated samples, harmonic values for each sample/node, ...)
    // the printed samples are parsed in python and used to create networkit graphs; in this way, mc is tested for correctness
    const uint64_t nodeCount = G.getNodeCount();
    const uint64_t edgeCount = G.getEdgeCount();
    std::vector<std::vector<double>> accumulatedResultsPerThread(OMP_MAX_THREADS, std::vector<double>(nodeCount, 0.0));

    std::vector<std::unordered_set<nodePair, nodePairHash>> outputEdges;
    outputEdges.resize(mcSampleCount);
    for (auto& set : outputEdges){
        set.reserve(edgeCount);
    }

    std::vector<std::vector<double>> singleNodeResultsPerRun(mcSampleCount);
    for (auto& run : singleNodeResultsPerRun){
        run.resize(nodeCount);
    }

    std::vector<std::unordered_map<nodePair, std::pair<double, double>, nodePairHash>> randNumVsProb(mcSampleCount);
    for (auto & map : randNumVsProb){
        map.reserve(edgeCount);
    }

#pragma omp parallel for num_threads(accumulatedResultsPerThread.size())
    for (uint64_t runNum = 0; runNum < mcSampleCount; ++runNum) {
        const auto threadId = static_cast<std::size_t>(omp_get_thread_num());
        std::random_device rand;
        std::default_random_engine eng(rand());
        std::uniform_real_distribution<double> distribution(0.0, 1.0);

        std::unordered_set<nodePair, nodePairHash> sampledEdges;
        sampledEdges.reserve(edgeCount);

        for (node current = 0; current < nodeCount - 1; ++current) {
            for (const node neighbor: G.neighbors(current)) {
                if (current > neighbor) {
                    continue;
                }
                const double probVal = G.getProb(current, neighbor);
                const double randVal = distribution(eng);
                if (randVal < probVal) {
                    sampledEdges.insert(nodePair{current, neighbor});
                    outputEdges[runNum].insert(nodePair{current, neighbor});
                }
                randNumVsProb[runNum].insert({{current, neighbor}, {randVal, probVal}});
            }
        }
        for (node v = 0; v < nodeCount; ++v){
            double harmonic = 0.0;
            std::vector<uint64_t> dist = allDistancesSingleSource(G, v, sampledEdges);
            for (node s = 0; s < v; ++s){
                if (dist[s] != INFINITY_UINT_64){
                    harmonic += 1.0 / static_cast<double>(dist[s]);
                }
            }
            for (node s = v+1; s < nodeCount; ++s){
                if (dist[s] != INFINITY_UINT_64){
                    harmonic += 1.0 / static_cast<double>(dist[s]);
                }
            }
            accumulatedResultsPerThread[threadId][v] += harmonic;
            singleNodeResultsPerRun[runNum][v] = harmonic / static_cast<double>(nodeCount-1);
        }
    }
    auto it = std::next(accumulatedResultsPerThread.begin());
    std::vector<double> mcResult(std::move(accumulatedResultsPerThread[0]));
    while (it != accumulatedResultsPerThread.end()){
        std::transform(it -> begin(), it -> end(),
                       mcResult.begin(), mcResult.begin(),
                       std::plus<>());
        ++it;
    }
    const double normalizationAndAveragingFactor = 1.0 / (static_cast<double>((nodeCount - 1)*mcSampleCount));
    std::transform(mcResult.begin(), mcResult.end(), mcResult.begin(),
                   [normalizationAndAveragingFactor](double elem){return elem * normalizationAndAveragingFactor;});
    for (uint64_t runNum = 0; runNum < mcSampleCount; ++runNum) {
        mcSampleString += "sample " + std::to_string(runNum + 1) + ":\n";
        for (node v = 0; v < nodeCount - 1; ++v) {
            for (node u = v + 1; u < nodeCount; ++u) {
                if (!G.hasEdge(u, v)) {
                    continue;
                }
                auto[rand, prob] = randNumVsProb[runNum].at({v, u});
                if (outputEdges[runNum].find({v, u}) != outputEdges[runNum].end()) {
                    mcSampleString += "\t(" + std::to_string(v) + ", " + std::to_string(u) + ") is sampled...";
                } else {
                    mcSampleString += "\t(" + std::to_string(v) + ", " + std::to_string(u) + ") is NOT sampled...";
                }
                mcSampleString += "rand: " + std::to_string(rand) + ", prob: " + std::to_string(prob) + "\n";
            }
        }
        for (node v = 0; v < nodeCount; ++v) {
            mcSampleString +=
                    "\t\tH(" + std::to_string(v) + ")= " + std::to_string(singleNodeResultsPerRun[runNum][v]) + "\n";
        }
    }
    return mcResult;
}

std::pair<std::vector<std::pair<double, std::vector<bool>>>, std::vector<nodePair>> allSpIgnoringUncertaintyBetweenness(
        const UncertainGraph &G, const node source, const node target, const std::unordered_set<nodePair, nodePairHash> &edgeRemoved){
    const uint64_t nodeCount = G.getNodeCount();

    std::vector<uint64_t> dist(nodeCount, INFINITY_UINT_64);
    dist[source] = 0;

    // pathsAndProbabilities[v] = vector of pairs (Pr(S), path_S) for every shortest path S from source to v
    // for S and every node v: path_S[v] = true if v was traversed on S
    std::vector<std::vector<std::pair<double, std::vector<bool>>>> pathsAndProbabilities(nodeCount);
    pathsAndProbabilities[source].emplace_back(std::make_pair(1.0, std::vector<bool>(nodeCount, false)));

    // minEdges[v] = (edge e, probability p, depth d)
    // at any point, e is an edge with minimal prob. p on all shortest paths source -> v, that was found most recently (i.e. at max. depth d)
    std::vector<std::tuple<nodePair, double, uint64_t>> minEdges(nodeCount, std::make_tuple(nodePair{0, 0}, INFINITY_DOUBLE, 0));

    std::queue<node> Q{};
    Q.push(source);
    while(!Q.empty()){
        const node current = Q.front();
        if (current == target){
            break;
        }
        Q.pop();
        for (const node child : G.neighbors(current)){
            if (edgeRemoved.find(child < current ? nodePair{child, current} : nodePair{current, child}) != edgeRemoved.end()){
                continue;
            }
            if (dist[child] == INFINITY_UINT_64){
                Q.push(child);
                dist[child] = dist[current] + 1;

                // check if {current,child} is minimal prob edge for all (current) shortest source -> child paths
                // else minEdges[child] = minEdges[current]
                double newProb = G.getProb(current, child);
                const double prevMinProb = std::get<1>(minEdges[current]);
                minEdges[child] = prevMinProb >= newProb ?
                        std::make_tuple(nodePair{current, child}, newProb, dist[child])
                        :  minEdges[current];

                // child inherits all shortest paths from current (i.e. their existence probabilities and inner nodes)
                // their probabilities are multiplied with P(current,child) though and
                // additionally, path[child]=true is set, as the given path has just traversed child
                pathsAndProbabilities[child].reserve(pathsAndProbabilities[current].size());
                std::transform(
                        pathsAndProbabilities[current].begin(),
                        pathsAndProbabilities[current].end(),
                        std::back_inserter(pathsAndProbabilities[child]),
                        [newProb, child](const std::pair<double, std::vector<bool>> &elem) -> std::pair<double, std::vector<bool>>
                        { std::vector<bool> newPath(elem.second);
                            newPath[child] = true;
                            return std::make_pair(elem.first * newProb, std::move(newPath));}
                );
            }
            else if (dist[child] == dist[current] + 1){
                // compare probabilities and (if needed) depth to determine min edge
                double newProb = G.getProb(current, child);
                const double prevMinProbChild = std::get<1>(minEdges[child]);
                const double prevMinProbCurrent = std::get<1>(minEdges[current]);

                if (std::min(prevMinProbChild, prevMinProbCurrent) >= newProb) {
                    minEdges[child] = std::make_tuple(nodePair{current, child}, newProb, dist[child]);
                } else if (prevMinProbChild > prevMinProbCurrent) {
                    minEdges[child] = minEdges[current];
                } else if (prevMinProbChild == prevMinProbCurrent) {
                    // if dist(source -> min_edge_current) > dist(source -> min_edge_child)
                    if (std::get<2>(minEdges[current]) > std::get<2>(minEdges[child])) {
                        minEdges[child] = minEdges[current];
                    }
                }

                // increase size of pathsAndProbabilities[child] to insert all paths that lead over current
                pathsAndProbabilities[child].reserve(pathsAndProbabilities[child].size() + pathsAndProbabilities[current].size());
                // insert the paths leading over current, but multiply the prob with P(current,child) and set path[child]=true
                std::transform(
                        pathsAndProbabilities[current].begin(),
                        pathsAndProbabilities[current].end(),
                        std::back_inserter(pathsAndProbabilities[child]),
                        [newProb, child](const std::pair<double, std::vector<bool>> &elem) -> std::pair<double, std::vector<bool>>
                        { std::vector<bool> newPath = elem.second;
                            newPath[child] = true;
                            return std::make_pair(elem.first * newProb, newPath);}
                );
            }
        }
    }
    // source, target disconnected; return empty result
    if (dist[target] == INFINITY_UINT_64){
        return std::make_pair(std::vector<std::pair<double, std::vector<bool>>>{}, std::vector<nodePair>{});
    }
    // edge {source,target} exists; no need to retrieve E_min ( = {{source,target}})
    if (dist[target] == 1){
        return std::make_pair(std::move(pathsAndProbabilities[target]), std::vector<nodePair>{nodePair{source, target}});
    }

    std::vector<nodePair> edgesToRemove; // no pre allocation for E_min. generally we will delete much fewer edges than paths in big graphs.
    std::vector<bool> visited(nodeCount, false);
    visited[target] = true;

    // delete {target, child} if the edge is minimal (and still exists)
    // else enqueue child, if edge {target,child} is not removed and child is predecessor on shortest path(s)
    Q = std::queue<node>{};
    for (const node child: G.neighbors(target)) {
        if (edgeRemoved.find(target < child ? nodePair{target, child} : nodePair{child, target}) != edgeRemoved.end()){
            continue;
        }
        if (dist[child] == dist[target] - 1) {
            visited[child] = true;
            if (G.getProb(target, child) <= std::get<1>(minEdges[child])) {
                edgesToRemove.emplace_back(target, child);
            } else {
                Q.push(child);
            }
        }
    }
    // traverse shortest paths backwards until min edge is found again; then stop traversal there and add edge to E_min
    while(!Q.empty()){
        const node current = Q.front();
        Q.pop();

        for (const node child : G.neighbors(current)){
            if (edgeRemoved.find(child < current ? nodePair{child, current} : nodePair{current, child}) != edgeRemoved.end()){
                continue;
            }
            if (dist[child] == dist[current] - 1 ){
                const auto [u, v] = std::get<0>(minEdges[current]);
                if (u == child && v == current){
                    edgesToRemove.emplace_back(u, v);
                }
                else if (!visited[child]){
                    visited[child] = true;
                    Q.push(child);
                }
            }
        }
    }
    return std::make_pair(std::move(pathsAndProbabilities[target]), std::move(edgesToRemove));
}

std::pair<std::vector<std::pair<double, std::vector<node>>>, std::vector<nodePair>> allSpIgnoringUncertaintyBetweennessMemory(
        const UncertainGraph &G, const node source, const node target, const std::unordered_set<nodePair, nodePairHash> &edgeRemoved){
    // i only comment on the differences to allSpIgnoringUncertaintyBetweenness(...), more detail is there
    const uint64_t nodeCount = G.getNodeCount();
    std::vector<uint64_t> dist(nodeCount, INFINITY_UINT_64);
    dist[source] = 0;

    // current maximal depth of dequeued nodes, init. zero
    uint64_t currentMaxDepth = 0;
    // vector holds all nodes that have been found at currentMaxDepth (init. just source)
    std::vector<node> nodesWithCurrentDepth;
    nodesWithCurrentDepth.push_back(source);

    std::vector<std::vector<std::pair<double, std::vector<node>>>> pathsAndProbabilities(nodeCount);
    pathsAndProbabilities[source].emplace_back(std::make_pair(1.0, std::vector<node>()));

    std::vector<std::tuple<nodePair, double, uint64_t>> minEdges(nodeCount, std::make_tuple(nodePair{0, 0}, INFINITY_DOUBLE, 0));

    std::queue<node> Q{};
    Q.push(source);
    while(!Q.empty()){
        const node current = Q.front();
        if (current == target){
            break;
        }
        Q.pop();
        // node with depth k+1 is dequeued for the first time; release the paths and probabilities from all nodes with depth k
        if (dist[current] > currentMaxDepth){
            for (node v : nodesWithCurrentDepth){
                pathsAndProbabilities[v].clear();
                pathsAndProbabilities[v].shrink_to_fit();
            }
            // set new currentMaxDepth and vector nodeWithCurrentDepth
            currentMaxDepth = dist[current];
            nodesWithCurrentDepth.clear();
            nodesWithCurrentDepth.push_back(current);
        }
        for (const node child : G.neighbors(current)){
            if (edgeRemoved.find(child < current ? nodePair{child, current} : nodePair{current, child}) != edgeRemoved.end()){
                continue;
            }
            if (dist[child] == INFINITY_UINT_64){
                Q.push(child);
                dist[child] = dist[current] + 1;
                double newProb = G.getProb(current, child);
                const double prevMinProb = std::get<1>(minEdges[current]);
                minEdges[child] = prevMinProb >= newProb ?
                        std::make_tuple(nodePair{current, child}, newProb, dist[child])
                        :  minEdges[current];

                pathsAndProbabilities[child].reserve(pathsAndProbabilities[current].size());
                std::transform(
                        pathsAndProbabilities[current].begin(),
                        pathsAndProbabilities[current].end(),
                        std::back_inserter(pathsAndProbabilities[child]),
                        [newProb, child](const std::pair<double, std::vector<node>> &elem) -> std::pair<double, std::vector<node>>
                        { std::vector<node> newPath(elem.second);
                            newPath.push_back(child);
                            return std::make_pair(elem.first * newProb, std::move(newPath));}
                );
            }
            else if (dist[child] == dist[current] + 1){
                double newProb = G.getProb(current, child);
                const double prevMinProbChild = std::get<1>(minEdges[child]);
                const double prevMinProbCurrent = std::get<1>(minEdges[current]);
                pathsAndProbabilities[child].reserve(pathsAndProbabilities[child].size() + pathsAndProbabilities[current].size());
                std::transform(
                        pathsAndProbabilities[current].begin(),
                        pathsAndProbabilities[current].end(),
                        std::back_inserter(pathsAndProbabilities[child]),
                        [newProb, child](const std::pair<double, std::vector<node>> &elem) -> std::pair<double, std::vector<node>>
                        { std::vector<node> newPath = elem.second;
                            newPath.push_back(child);
                            return std::make_pair(elem.first * newProb, newPath);}
                );
                if (std::min(prevMinProbChild, prevMinProbCurrent) >= newProb) {
                    minEdges[child] = std::make_tuple(nodePair{current, child}, newProb, dist[child]);
                } else if (prevMinProbChild > prevMinProbCurrent) {
                    minEdges[child] = minEdges[current];
                } else if (prevMinProbChild == prevMinProbCurrent) {
                    // if dist(source -> min_edge_current) > dist(source -> min_edge_child)
                    if (std::get<2>(minEdges[current]) > std::get<2>(minEdges[child])) {
                        minEdges[child] = minEdges[current];
                    }
                }
            }
        }
    }
    if (dist[target] == INFINITY_UINT_64){
        return std::make_pair(std::vector<std::pair<double, std::vector<node>>>{}, std::vector<nodePair>{});
    }
    if (dist[target] == 1){
        return std::make_pair(std::move(pathsAndProbabilities[target]), std::vector<nodePair>{nodePair{source, target}});
    }
    std::vector<nodePair> edgesToRemove;
    std::vector<bool> visited(nodeCount, false);
    visited[target] = true;

    Q = std::queue<node>{};
    for (const node child: G.neighbors(target)) {
        if (edgeRemoved.find(target < child ? nodePair{target, child} : nodePair{child, target}) != edgeRemoved.end()){
            continue;
        }
        if (dist[child] == dist[target] - 1) {
            visited[child] = true;
            if (G.getProb(target, child) <= std::get<1>(minEdges[child])) {
                edgesToRemove.emplace_back(target, child);
            } else {
                Q.push(child);
            }
        }
    }
    while(!Q.empty()){
        const node current = Q.front();
        Q.pop();

        for (const node child : G.neighbors(current)){
            if (edgeRemoved.find(child < current ? nodePair{child, current} : nodePair{current, child}) != edgeRemoved.end()){
                continue;
            }
            if (dist[child] == dist[current] - 1 ){
                const auto [u, v] = std::get<0>(minEdges[current]);
                if (u == child && v == current){
                    edgesToRemove.emplace_back(u, v);
                }
                else if (!visited[child]){
                    visited[child] = true;
                    Q.push(child);
                }
            }
        }
    }
    return std::make_pair(std::move(pathsAndProbabilities[target]), std::move(edgesToRemove));
}

std::vector<double> parallelBetweennessPSP(const UncertainGraph &G, const double phi) {
    const uint64_t nodeCount = G.getNodeCount();

    // init. partial betweenness sums of |V| zeros per thread
    std::vector<std::vector<double>> betweennessResultPerThread(OMP_MAX_THREADS, std::vector<double>(nodeCount, 0.0));

    // loop over s,t in parallel to calculate PSP betweenness for each such node pair
#pragma omp parallel for num_threads(betweennessResultPerThread.size())
    for (node s = 0; s < nodeCount - 1; ++s){
        std::vector<double> currentThreadResult(nodeCount, 0.0); // partial PSP betweenness sum for nodes s,t
        for (node t = s+1; t < nodeCount; ++t){
            // reset partial sum after t changes
            std::fill(currentThreadResult.begin(), currentThreadResult.end(), 0.0);
            const auto threadId = static_cast<std::size_t>(omp_get_thread_num());

            std::unordered_set<nodePair, nodePairHash> removedEdges;
            removedEdges.reserve(G.getEdgeCount());


            double relativeProbSum = 0.0;                   // sum of all est. rel. probabilities of s,t paths


            double negatedProbabilityProduct = 1.0;         // negatedProbabilityProduct = prod_{p in P_st} (1-Pr(p)) at any given point
            double phiST = 0.0;                             // est. connection probability phi_st

            while (phiST < phi){
                // explore next set of paths
                auto [pathsAndProbabilities, edgesToRemove] = allSpIgnoringUncertaintyBetweenness(G, s, t, removedEdges);
                if (pathsAndProbabilities.empty()){
                    // empty result -> s,t disconnected
                    break;
                }

                // we need the current negatedProbProd to calculate the est. rel. prob. of S = Pr(S)*negatedProbabilityProduct
                // yet we also want to calculate negatedProbabilityProduct *= (Pr(s)), hence the temporary variable
                double newNegatedProbabilityProduct = negatedProbabilityProduct;

                for (auto &[absolutePathProb, isOnPath] : pathsAndProbabilities){
                    // Pr(S)*negatedProbabilityProduct
                    const double newRelativePathProb = absolutePathProb * negatedProbabilityProduct;

                    // update relativeProbSum and (temporary variable) newNegatedProbabilityProduct
                    relativeProbSum += newRelativePathProb;
                    newNegatedProbabilityProduct *= (1.0 - absolutePathProb);

                    // add new est. relative path probability to all inner nodes on current path
                    for (node innerNode = 0; innerNode < s; ++innerNode){
                        if (isOnPath[innerNode]){
                            currentThreadResult[innerNode] += newRelativePathProb;
                        }
                    }
                    for (node innerNode = s+1; innerNode < t; ++innerNode){
                        if (isOnPath[innerNode]){
                            currentThreadResult[innerNode] += newRelativePathProb;
                        }
                    }
                    for (node innerNode = t+1; innerNode < nodeCount; ++innerNode){
                        if (isOnPath[innerNode]){
                            currentThreadResult[innerNode] += newRelativePathProb;
                        }
                    }
                }
                // discard of temporary newNegatedProbabilityProduct; update phiST; remove all edges from edgesToRemove
                negatedProbabilityProduct = newNegatedProbabilityProduct;
                phiST = 1.0 - negatedProbabilityProduct;
                for (auto &[u,v] : edgesToRemove){
                    removedEdges.insert(u < v ? nodePair{u, v} : nodePair{v, u});
                }
            }
            if (relativeProbSum == 0.0){
                // no s,t paths
                continue;
            }
            // every result of nodes v != s,t has to be multiplied by phiST/relativeProbSum according to PSP betweenness formula
            // the result for each node v != s,t is then added to the partial betweenness sum of the given thread
            const double phiSTDivByRelProbSum = phiST / static_cast<double>(relativeProbSum);
            for (node v = 0; v < s; ++v){
                betweennessResultPerThread[threadId][v] += currentThreadResult[v] * phiSTDivByRelProbSum;
            }
            for (node v = s+1; v < t; ++v){
                betweennessResultPerThread[threadId][v] += currentThreadResult[v] * phiSTDivByRelProbSum;
            }
            for (node v = t+1; v < nodeCount; ++v){
                betweennessResultPerThread[threadId][v] += currentThreadResult[v] * phiSTDivByRelProbSum;
            }
        }
    }
    // after iteration of all s,t nodes is finished, the results of all threads are accumulated
    // move first partial result vector, then add the remaining vectors
    auto it = std::next(betweennessResultPerThread.begin());
    std::vector<double> betweennessResult(std::move(betweennessResultPerThread[0]));
    while (it != betweennessResultPerThread.end()){
        std::transform(it -> begin(), it -> end(),
                       betweennessResult.begin(), betweennessResult.begin(), std::plus<>());
        ++it;
    }
    // finally, the accumulated results are normalized with 2 / ((n-1)(n-2))
    const double normalizationFactor = 2.0 / (static_cast<double>(nodeCount-1) * static_cast<double>(nodeCount-2));
    std::transform(betweennessResult.begin(), betweennessResult.end(), betweennessResult.begin(),
                   [normalizationFactor](const double elem) {return normalizationFactor * elem;});
    return betweennessResult;
}

std::vector<double> parallelBetweennessPSPMemory(const UncertainGraph &G, const double phi) {
    // the only difference to parallelBetweennessPSP(...) is that allSpIgnoringUncertaintyBetweennessMemory(...) gets called
    const uint64_t nodeCount = G.getNodeCount();
    std::vector<std::vector<double>> betweennessResultPerThread(OMP_MAX_THREADS, std::vector<double>(nodeCount, 0.0));

#pragma omp parallel for num_threads(betweennessResultPerThread.size())
    for (node s = 0; s < nodeCount - 1; ++s){
        std::vector<double> currentThreadResult(nodeCount, 0.0);
        for (node t = s+1; t < nodeCount; ++t){
            std::fill(currentThreadResult.begin(), currentThreadResult.end(), 0.0);
            const auto threadId = static_cast<std::size_t>(omp_get_thread_num());
            std::unordered_set<nodePair, nodePairHash> removedEdges;
            removedEdges.reserve(G.getEdgeCount());


            double relativeProbSum = 0.0;
            double negatedProbabilityProduct = 1.0;
            double phiST = 0.0;

            while (phiST < phi){
                auto [pathsAndProbabilities, edgesToRemove] =
                        allSpIgnoringUncertaintyBetweennessMemory(G, s, t, removedEdges);
                if (pathsAndProbabilities.empty()){
                    break;
                }
                double newNegatedProbabilityProduct = negatedProbabilityProduct;
                for (auto &[absolutePathProb, path] : pathsAndProbabilities){
                    const double newRelativePathProb = absolutePathProb * negatedProbabilityProduct;
                    relativeProbSum += newRelativePathProb;
                    newNegatedProbabilityProduct *= (1.0 - absolutePathProb);
                    for (node innerNode : path){
                        if (innerNode == t){
                            continue;
                        }
                        currentThreadResult[innerNode] += newRelativePathProb;
                    }
                }
                negatedProbabilityProduct = newNegatedProbabilityProduct;
                phiST = 1.0 - negatedProbabilityProduct;
                for (auto &[u,v] : edgesToRemove){
                    removedEdges.insert(u < v ? nodePair{u, v} : nodePair{v, u});
                }
            }
            if (relativeProbSum == 0.0){
                continue;
            }
            const double phiSTDivByRelProbSum = phiST / static_cast<double>(relativeProbSum);
            for (node v = 0; v < s; ++v){
                betweennessResultPerThread[threadId][v] += currentThreadResult[v] * phiSTDivByRelProbSum;
            }
            for (node v = s+1; v < t; ++v){
                betweennessResultPerThread[threadId][v] += currentThreadResult[v] * phiSTDivByRelProbSum;
            }
            for (node v = t+1; v < nodeCount; ++v){
                betweennessResultPerThread[threadId][v] += currentThreadResult[v] * phiSTDivByRelProbSum;
            }
        }
    }
    auto it = std::next(betweennessResultPerThread.begin());
    std::vector<double> betweennessResult(std::move(betweennessResultPerThread[0]));
    while (it != betweennessResultPerThread.end()){
        std::transform(it -> begin(), it -> end(),
                       betweennessResult.begin(), betweennessResult.begin(), std::plus<>());
        ++it;
    }
    const double normalizationFactor = 2.0 / (static_cast<double>(nodeCount-1) * static_cast<double>(nodeCount-2));
    std::transform(betweennessResult.begin(), betweennessResult.end(), betweennessResult.begin(),
                   [normalizationFactor](const double elem) {return normalizationFactor * elem;});
    return betweennessResult;
}

std::vector<double>& brandesBetweennessAddToAllNodes(const UncertainGraph &G, const std::unordered_set<nodePair, nodePairHash>& sampledEdges, std::vector<double>& result) {
    // straightforward implementation of brandes algorithm for betweenness centrality of all nodes in G (using only the edges in sampledEdges).
    // source: A Faster Algorithm for Betweenness Centrality, Ulrik Brandes
    //      http://snap.stanford.edu/class/cs224w-readings/brandes01centrality.pdf
    // the betweenness is added to the result vector though, to avoid the creation of a temporary vector of size |V| for every monte carlo iteration
    const uint64_t nodeCount = G.getNodeCount();
    std::vector<node> reverseRetrievalStack;
    reverseRetrievalStack.reserve(nodeCount);
    std::vector<std::vector<node>>predecessorList(nodeCount, std::vector<node>{});
    std::vector<node> sigma(nodeCount, 0);
    std::vector<node> distance(nodeCount, INFINITY_UINT_64);

    for (node s = 0; s < nodeCount; ++s){
        reverseRetrievalStack.clear();
        for (auto& list : predecessorList){
            list.clear();
        }
        std::fill(sigma.begin(), sigma.end(), 0);
        sigma[s] = 1;
        std::fill(distance.begin(), distance.end(), INFINITY_UINT_64);
        distance[s] = 0;
        std::queue<node> Q{};
        Q.push(s);
        while(!Q.empty()){
            const node v = Q.front();
            Q.pop();
            reverseRetrievalStack.push_back(v);
            for (const node w : G.neighbors(v)){
                if (sampledEdges.find(v < w ? nodePair{v,w} : nodePair{w,v}) == sampledEdges.end()){
                    continue;
                }
                if (distance[w] == INFINITY_UINT_64){
                    Q.push(w);
                    distance[w] = distance[v] + 1;
                }
                if (distance[w] == distance[v] + 1){
                    sigma[w] = sigma[w] + sigma[v];
                    predecessorList[w].push_back(v);
                }
            }
        }
        std::vector<double> delta(nodeCount, 0.0);
        while (!reverseRetrievalStack.empty()){
            const node w = reverseRetrievalStack.back();
            reverseRetrievalStack.pop_back();
            for (const node v : predecessorList[w]){
                delta[v] += (static_cast<double>(sigma[v]) / static_cast<double>(sigma[w])) * (1 + delta[w]);
            }
            if (w != s){
                result[w] += delta[w];
            }
        }
    }
    return result;
}

std::vector<double> parallelMcBetweenness(const UncertainGraph& G, const uint64_t mcSampleCount){
    const uint64_t nodeCount = G.getNodeCount();
    const uint64_t edgeCount = G.getEdgeCount();

    // init partial betweenness sums; |V| zeros for each thread
    std::vector<std::vector<double>> accumulatedResultsPerThread(OMP_MAX_THREADS, std::vector<double>(nodeCount, 0.0));

    // sample mcSampleCount graphs, distribute them over all available OMP threads
#pragma omp parallel for num_threads(accumulatedResultsPerThread.size())
    for (uint64_t runNum = 0; runNum < mcSampleCount; ++runNum) {
        const auto threadId = static_cast<std::size_t>(omp_get_thread_num());

        std::random_device rand;
        std::default_random_engine eng(rand());
        std::uniform_real_distribution<double> distribution(0.0, 1.0);

        std::unordered_set<nodePair, nodePairHash> sampledEdges;
        sampledEdges.reserve(edgeCount);

        for (node current = 0; current < nodeCount - 1; ++current) {
            for (const node neighbor: G.neighbors(current)) {
                if (current > neighbor) { // dont sample edges twice
                    continue;
                }
                if (distribution(eng) < G.getProb(current, neighbor)) {
                    sampledEdges.insert(nodePair{current, neighbor});
                }
            }
        }
        // add betweenness of all nodes in current instance to partial result of this thread (using augmented brandes algorithm)
        accumulatedResultsPerThread[threadId] = brandesBetweennessAddToAllNodes(G, sampledEdges, accumulatedResultsPerThread[threadId]);
    }
    // iterator starts at second element; first element (partial sum with thread id 0) is moved into mcResult.
    auto threadResultIterator = std::next(accumulatedResultsPerThread.begin());
    std::vector<double> mcResult(std::move(accumulatedResultsPerThread[0]));
    // add partial results of all other threads to mcResult
    while (threadResultIterator != accumulatedResultsPerThread.end()){
        std::transform(threadResultIterator -> begin(), threadResultIterator -> end(),
                       mcResult.begin(), mcResult.begin(),std::plus<>());
        ++threadResultIterator;
    }
    // normalize and get MC average, i.e. multiply by 1 / (mcSamplecount(n-1)(n-2))
    // 1 / ... instead of  2 / ... because brandes already produces the factor two for undirected graphs
    const double normalizationAndAveragingFactor = 1.0 / static_cast<double>(mcSampleCount*(nodeCount-1)*(nodeCount-2));
    std::transform(mcResult.begin(), mcResult.end(), mcResult.begin(),
                   [normalizationAndAveragingFactor](double elem){return elem * normalizationAndAveragingFactor;});
    return mcResult;
}

std::vector<double> parallelMcBetweennessPrintSamples(const UncertainGraph& G, uint64_t mcSampleCount, std::string& mcSampleString){
    // this is just for debugging/testing; same as parallelMcBetweenness(...) with additional output (generated samples, betweenness values for each sample/node, ...)
    // the printed samples are parsed in python and used to create networkit graphs; in this way, mc is tested for correctness
    const uint64_t nodeCount = G.getNodeCount();
    const uint64_t edgeCount = G.getEdgeCount();
    std::vector<std::vector<double>> accumulatedResultsPerThread(OMP_MAX_THREADS, std::vector<double>(nodeCount, 0.0));

    std::vector<std::unordered_set<nodePair, nodePairHash>> outputEdges(mcSampleCount);
    for (auto& set : outputEdges){
        set.reserve(edgeCount);
    }

    std::vector<std::vector<double>> singleNodeResultsPerRun(mcSampleCount, std::vector<double>(nodeCount, 0.0));

    std::vector<std::unordered_map<nodePair, std::pair<double, double>, nodePairHash>> randNumVsProb(mcSampleCount);
    for (auto & map : randNumVsProb){
        map.reserve(edgeCount);
    }

#pragma omp parallel for num_threads(OMP_MAX_THREADS)
    for (uint64_t runNum = 0; runNum < mcSampleCount; ++runNum) {
        const auto threadId = static_cast<std::size_t>(omp_get_thread_num());
        std::random_device rand;
        std::default_random_engine eng(rand());
        std::uniform_real_distribution<double> distribution(0.0, 1.0);

        std::unordered_set<nodePair, nodePairHash> sampledEdges;
        sampledEdges.reserve(edgeCount);

        for (node current = 0; current < nodeCount - 1; ++current) {
            for (const node neighbor: G.neighbors(current)) {
                if (current > neighbor) {
                    continue;
                }
                const double probVal = G.getProb(current, neighbor);
                const double randVal = distribution(eng);
                if (randVal < probVal) {
                    sampledEdges.insert(nodePair{current, neighbor});
                    outputEdges[runNum].insert(nodePair{current, neighbor});
                }
                randNumVsProb[runNum].insert({{current, neighbor}, {randVal, probVal}});
            }
        }
        singleNodeResultsPerRun[runNum] = brandesBetweennessAddToAllNodes(G, sampledEdges, singleNodeResultsPerRun[runNum]);
        accumulatedResultsPerThread[threadId] = brandesBetweennessAddToAllNodes(G, sampledEdges, accumulatedResultsPerThread[threadId]);
    }
    auto threadResultIterator = std::next(accumulatedResultsPerThread.begin());
    std::vector<double> mcResult(std::move(accumulatedResultsPerThread[0]));
    while (threadResultIterator != accumulatedResultsPerThread.end()){
        std::transform(threadResultIterator -> begin(), threadResultIterator -> end(),
                       mcResult.begin(), mcResult.begin(), std::plus<>());
        ++threadResultIterator;
    }
    const double normalizationAndAveragingFactor = 1.0 / static_cast<double>(mcSampleCount*(nodeCount-1)*(nodeCount-2));
    std::transform(mcResult.begin(), mcResult.end(), mcResult.begin(),
                   [normalizationAndAveragingFactor](double elem){return elem * normalizationAndAveragingFactor;});

    for (uint64_t runNum = 0; runNum < mcSampleCount; ++runNum) {
        mcSampleString += "sample " + std::to_string(runNum + 1) + ":\n";
        for (node v = 0; v < nodeCount - 1; ++v) {
            for (node u = v + 1; u < nodeCount; ++u) {
                if (!G.hasEdge(u, v)) {
                    continue;
                }
                auto[rand, prob] = randNumVsProb[runNum].at({v, u});
                if (outputEdges[runNum].find({v, u}) != outputEdges[runNum].end()) {
                    mcSampleString += "\t(" + std::to_string(v) + ", " + std::to_string(u) + ") is sampled...";
                } else {
                    mcSampleString += "\t(" + std::to_string(v) + ", " + std::to_string(u) + ") is NOT sampled...";
                }
                mcSampleString += "rand: " + std::to_string(rand) + ", prob: " + std::to_string(prob) + "\n";
            }
        }
        const double normalizationFactor = 1.0 / static_cast<double>((nodeCount-1)*(nodeCount-2));
        for (node v = 0; v < nodeCount; ++v) {
            mcSampleString += "\t\tCB(" + std::to_string(v) + ")= " + std::to_string(singleNodeResultsPerRun[runNum][v] * normalizationFactor) + "\n";
        }
    }
    return mcResult;
}

std::string distanceTest(UncertainGraph& G, const double phi, const uint32_t OUTPUT_PRECISION){
    // just for testing/debugging purposes
    std::stringstream ss;
    std::string output = "psp distance distribution:\n";
    std::vector<double> distribution{};
    for (node s = 0; s < G.getNodeCount() - 1; ++s){
        for (node t = s + 1; t < G.getNodeCount(); ++t){
            distribution = pspDistanceDistribution(G, s, t, phi);
            for ( uint64_t i = 0 ; i < distribution.size() - 1; ++i){
                ss = std::stringstream();
                ss.precision(OUTPUT_PRECISION);
                ss << std::fixed << distribution[i];
                output += "\t\tp_" + std::to_string(s) + "," + std::to_string(t) + " (" + std::to_string(i+1) + ") = " + ss.str() + "\n";
            }
            ss = std::stringstream();
            ss.precision(OUTPUT_PRECISION);
            ss << std::fixed << distribution[distribution.size()-1];
            output += "\t\tp_" + std::to_string(s) + "," + std::to_string(t) + " (inf) = " + ss.str() + "\n";

            double dist;
            output += "\tdistance:\n";
            output += "\td_ER(" + std::to_string(s) + "," + std::to_string(t) + ") = ";
            ss = std::stringstream();
            ss.precision(OUTPUT_PRECISION);
            dist = getDistanceER(distribution);
            if (dist == INFINITY_DOUBLE){
                output += "inf\n";
            }
            else{
                ss << std::fixed << dist;
                output += ss.str() + "\n";
            }
            output += "\td_MED(" + std::to_string(s) + "," + std::to_string(t) + ") = ";
            ss = std::stringstream();
            ss.precision(OUTPUT_PRECISION);
            dist = getDistanceMed(distribution);
            if (dist == INFINITY_DOUBLE){
                output += "inf\n";
            }
            else{
                ss << std::fixed << dist;
                output += ss.str() + "\n";
            }
            output += "\td_MAJ(" + std::to_string(s) + "," + std::to_string(t) + ") = ";
            ss = std::stringstream();
            ss.precision(OUTPUT_PRECISION);
            dist = getDistanceMaj(distribution);
            if (dist == INFINITY_DOUBLE){
                output += "inf\n";
            }
            else{
                ss << std::fixed << dist;
                output += ss.str() + "\n";
            }
        }
    }
    return output;
}
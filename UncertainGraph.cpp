#include "UncertainGraph.hpp"
#include <algorithm>
#include <cassert>

UncertainGraph::UncertainGraph(uint64_t n, std::vector<edge>& edgeList) : UncertainGraph(n, std::as_const(edgeList)){}

UncertainGraph::UncertainGraph(const uint64_t n, const std::vector<edge> &edgeList)
        : nodeCount(n), edgeCount(edgeList.size()), E(n) {
    P.reserve(edgeList.size());
    for (const edge & e : edgeList){
        assert((e.u < n) && (e.v < n));
        if (hasEdge(e.u, e.v) || e.probability == 0.0){
            --edgeCount;
            continue;
        }
        // we need both vectors E[u] and E[v] to get fast access to neighbors of any node
        E[e.u].push_back(e.v);
        E[e.v].push_back(e.u);
        P.insert({e.u < e.v ? nodePair{e.u, e.v} : nodePair{e.v, e.u}, e.probability});
    }
    for (std::vector<node>& list : E){
        list.shrink_to_fit();
    }
}

// copy constructor, copy assignment, move constructor, move assignment and destructor (=default in hpp)
// all are trivial implementations and generally not needed (just here for best practice; of rule of five)
UncertainGraph::UncertainGraph(const UncertainGraph &other) {
    nodeCount = other.nodeCount;
    edgeCount = other.edgeCount;
    E = other.E;
    P = other.P;
}

UncertainGraph& UncertainGraph::operator=(const UncertainGraph &other) {
    nodeCount = other.nodeCount;
    edgeCount = other.edgeCount;
    E = other.E;
    P = other.P;
    return *this;
}

UncertainGraph::UncertainGraph(UncertainGraph &&other) noexcept {
    nodeCount = other.nodeCount;
    edgeCount = other.edgeCount;
    E = std::move(other.E);
    P = std::move(other.P);
}

UncertainGraph& UncertainGraph::operator=(UncertainGraph &&other) noexcept {
    nodeCount = other.nodeCount;
    edgeCount = other.edgeCount;
    E = std::move(other.E);
    P = std::move(other.P);
    return *this;
}

bool UncertainGraph::hasEdge(const node u, const node v) const noexcept{
    auto it = std::find(E[u].begin(), E[u].end(), v);
    return it != E[u].end();
}

double UncertainGraph::getProb(const node u, const node v) const{
    return P.at(u < v ? nodePair{u,v} : nodePair{v,u});
}

uint64_t UncertainGraph::getNodeCount() const noexcept{
    return nodeCount;
}

uint64_t UncertainGraph::getEdgeCount() const noexcept{
    return edgeCount;
}

const std::vector<node>& UncertainGraph::neighbors(node v) const noexcept{
    return E[v];
}
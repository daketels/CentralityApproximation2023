#ifndef UNCERTAINCENTRALITY_UNCERTAINGRAPH_HPP
#define UNCERTAINCENTRALITY_UNCERTAINGRAPH_HPP

#include <unordered_map>
#include <vector>
#include <utility>
#include <boost/functional/hash.hpp>

using node = uint64_t;
using nodePair = std::pair<node, node>;
using nodePairHash = boost::hash<nodePair>;

struct edge {
    node u;
    node v;
    double probability;
    edge(node u, node v, double probability) : u(u), v(v), probability(probability) {};
    edge() = default;
};

class UncertainGraph {
private:
    // a graph with nodeCount = n includes the nodes 0, ... , n-1
    uint64_t nodeCount;
    uint64_t edgeCount;

    // for each edge {u, v, p}
    //      v is saved in E[u], u is saved in E[v] (creating both entries greatly improves performance when retrieving all neighbors of a node)
    //      p is saved in P[min(u,v), max(u,v)]
    std::vector<std::vector<node>> E;
    std::unordered_map<nodePair, double, nodePairHash> P{};
public:
    UncertainGraph() = delete;

    // Constructor delegates to UncertainGraph(uint64_t n, const std::vector<edge>& edgeList)
    explicit UncertainGraph(uint64_t n, std::vector<edge>& edgeList);

    // construct graph with n nodes and all edges in edgeList
    // edges with probability 0 and loops are ignored
    explicit UncertainGraph(uint64_t n, const std::vector<edge>& edgeList);

    UncertainGraph(const UncertainGraph& other);

    UncertainGraph& operator=(const UncertainGraph& other);

    UncertainGraph(UncertainGraph&& other) noexcept;

    UncertainGraph& operator=(UncertainGraph&& other) noexcept;

    ~UncertainGraph() = default;

    // hasEdge(u, v) returns true iff there is an edge between u and v
    [[nodiscard]] bool hasEdge(node u, node v) const noexcept;

    // returns probability of the edge {u, v}
    [[nodiscard]] double getProb(node u, node v) const;

    // returns number of nodes in the graph
    [[nodiscard]] uint64_t getNodeCount() const noexcept;

    // returns number of edges in the graph
    [[nodiscard]] uint64_t getEdgeCount() const noexcept;

    // returns const ref. to vector of all neighbors of node v
    [[nodiscard]] const std::vector<node>& neighbors(node v) const noexcept;
};
#endif //UNCERTAINCENTRALITY_UNCERTAINGRAPH_HPP
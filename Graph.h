// Graph.h
#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <queue>
#include <utility>
#include <limits>

class Graph {
public:
    Graph(int n);
    void addEdge(int u, int v, double cost);
    void initialize(int n);
    void dijkstra(int src, std::vector<double>& distMin);

private:
    struct Edge {
        int v;
        double cost;
    };
    std::vector<std::vector<Edge>> adjList;
    int numNodes;
};

#endif // GRAPH_H

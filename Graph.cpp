// Graph.cpp
#include "Graph.h"

Graph::Graph(int n) {
    initialize(n);
}

void Graph::initialize(int n) {
    numNodes = n;
    adjList.resize(n);
}

void Graph::addEdge(int u, int v, double cost) {
    std::cout << "Adding edge from node " << u << " to node " << v << std::endl;
    adjList[u].push_back({v, cost});
    adjList[v].push_back({u, cost}); // Assuming it's an undirected graph
}

void Graph::dijkstra(int src, std::vector<double>& distMin) {
    distMin.assign(numNodes, std::numeric_limits<double>::infinity());
    distMin[src] = 0;
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> pq;
    pq.push({0, src});

    while (!pq.empty()) {
        int u = pq.top().second;
        double distU = pq.top().first;
        pq.pop();

        if (distU > distMin[u]) continue;

        for (const auto& edge : adjList[u]) {
            int v = edge.v;
            double cost = edge.cost;
            if (distMin[u] + cost < distMin[v]) {
                distMin[v] = distMin[u] + cost;
                pq.push({distMin[v], v});
            }
        }
    }
}

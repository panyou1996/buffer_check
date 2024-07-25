#include <queue>
#include <vector>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <string>

#define FIN "dijkstra_test.in"
#define FOUT "dijkstra.out"
#define INF ((1LL << 31) - 1)
#define MAXN 500000

using namespace std;

struct Node {
    int y, cost;
    Node* next;
};

Node* V[MAXN];
bool inQueue[MAXN];
priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
int distMin[MAXN];
int nodes, edges;

void addEdge(const int x, const int y, const int cost) {
    Node* node = new Node;
    node->y = y;
    node->cost = cost;
    node->next = V[x];
    V[x] = node;
}

void readInput() {
    int x, y, c;
    FILE* fin = fopen(FIN, "r");

    if (fin == NULL) {
        perror("Error opening input file");
        exit(EXIT_FAILURE);
    }

    fscanf(fin, "%d %d", &nodes, &edges);

    for (int i = 0; i < edges; ++i) {
        fscanf(fin, "%d %d %d", &x, &y, &c);
        addEdge(x, y, c);
    }

    fclose(fin);
}

void initialize(int n) {
    for (int i = 0; i < n; ++i) {
        V[i] = nullptr;
        inQueue[i] = false;
        distMin[i] = INF;
    }
}

void dijkstra(int src) {
    distMin[src] = 0;
    pq.push({0, src});

    while (!pq.empty()) {
        int x = pq.top().second;
        pq.pop();

        if (inQueue[x]) continue;
        inQueue[x] = true;

        for (Node* p = V[x]; p != nullptr; p = p->next) {
            if (distMin[p->y] > distMin[x] + p->cost) {
                distMin[p->y] = distMin[x] + p->cost;
                pq.push({distMin[p->y], p->y});
            }
        }
    }
}

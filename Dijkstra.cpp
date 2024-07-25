#include <queue>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <string>
// 定义输入输出文件名
#define FIN "dijkstra_test.in"
#define FOUT "dijkstra.out"
// 定义无穷大的值
#define INF ((1LL << 31) - 1)
// 定义最大节点数
#define MAXN 500000

using namespace std;

// 定义图的节点结构
struct Node {
    int y, cost; // 目标节点和边的权重
    Node *next;   // 下一个邻接节点的指针
};

// 邻接表数组，存储每个节点的邻接节点
Node *V[MAXN];
// 标记节点是否在队列中
bool inQueue[MAXN];
// 优先队列，存储待处理的节点
priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
// 存储从起点到每个节点的最短距离
int distMin[MAXN];

// 定义节点数和边数
int nodes, edges;

// 添加边的函数
void addEdge(const int x, const int y, const int cost) {
    Node *node = new Node;
    node->y = y;
    node->cost = cost;
    node->next = V[x];
    V[x] = node;
}

// 读取输入文件的函数
void readInput() {
    int x, y, c;
    FILE *fin = fopen(FIN, "r");

    if (fin == NULL) {
        perror("Error opening input file");
        exit(EXIT_FAILURE);
    }

    // 读取节点数和边数
    fscanf(fin, "%d %d", &nodes, &edges);

    // 读取每条边的信息，并添加到图中
    for (int i = 0; i < edges; ++i) {
        fscanf(fin, "%d %d %d", &x, &y, &c);
        addEdge(x, y, c);
    }

    fclose(fin);
}

// Dijkstra算法求解最短路径的函数
void dijkstra(int start, int end) {
    // 初始化最短距离数组
    for (int i = 1; i <= nodes; i++) distMin[i] = INF;
    distMin[start] = 0;

    // 将起点加入优先队列
    pq.push({0, start});
    inQueue[start] = true;

    // 循环直到队列为空
    while (!pq.empty()) {
        auto [currentCost, currentNode] = pq.top();
        pq.pop();
        inQueue[currentNode] = false;

        // 遍历当前节点的所有邻接节点
        for (Node *p = V[currentNode]; p; p = p->next) {
            // 如果找到了更短的路径，则更新距离并加入队列
            if (distMin[p->y] > currentCost + p->cost) {
                distMin[p->y] = currentCost + p->cost;
                if (!inQueue[p->y]) {
                    pq.push({distMin[p->y], p->y});
                    inQueue[p->y] = true;
                }
            }
        }
    }

    // 输出从起始节点到目标节点的最短距离
    printf("The shortest distance from node %d to node %d is %s\n", start, end, distMin[end] < INF ? to_string(distMin[end]).c_str() : "Infinity");
}

// 写入输出结果的函数
void writeOutput(int start, int end) {
    FILE *fout = fopen(FOUT, "w");

    if (fout == NULL) {
        perror("Error opening output file");
        exit(EXIT_FAILURE);
    }

    dijkstra(start, end);

    // 写入从起始节点到目标节点的最短距离
    fprintf(fout, "%d\n", distMin[end] < INF ? distMin[end] : 0);

    // 关闭输出文件
    fclose(fout);
}

// 主函数
int main() {
    int startNode, endNode;
    startNode = 1;
    endNode = 4;
    readInput();
    writeOutput(startNode, endNode);
    return 0;
}
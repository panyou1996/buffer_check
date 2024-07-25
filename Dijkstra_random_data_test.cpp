#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>

#define EDGES 10000
#define MAX_VERTICES 5000

using namespace std;

int main() {
    // 初始化随机数生成器
    srand(static_cast<unsigned>(time(0)));

    // 创建输出文件流
    //ofstream outfile("dijkstra_test.in");
    ofstream outfile("input.txt");
    // 随机生成顶点数，这里假设最大顶点数为5000
    int vertices = 5000;
    outfile << vertices << " " << EDGES << 1 << endl;

    // 生成边
    for (int i = 0; i < EDGES; ++i) {
        int from = rand() % vertices + 1; // 顶点编号从1开始
        int to = rand() % vertices + 1;   // 顶点编号从1开始
        int weight = rand() % 100 + 1;   // 权重在1到100之间

        outfile << from << " " << to << " " << weight << endl;
    }

    // 关闭文件流
    outfile.close();

    cout << "Test data file 'dijkstra_test.in' generated successfully." << endl;

    return 0;
}

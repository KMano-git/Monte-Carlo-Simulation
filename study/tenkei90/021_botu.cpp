#include <bits/stdc++.h>
using namespace std;

// memo[start][goal]: -1=未計算, 0=到達不可, 1=到達可能

bool dfs(int cur, int goal, vector<vector<int>> &graph, vector<bool> &visited) {
    if (cur == goal) return true;
    if (visited[cur]) return false;
    visited[cur] = true;
    
    for (int next : graph[cur]) {
        if (dfs(next, goal, graph, visited)) return true;
    }
    return false;
}

// メモ化付きのcanreach（呼び出し側でキャッシュ）
bool canreach(int start, int goal, vector<vector<int>> &graph) {
    if (memo[start][goal] != -1) return memo[start][goal];
    
    vector<bool> visited(graph.size(), false);
    bool result = dfs(start, goal, graph, visited);
    memo[start][goal] = result ? 1 : 0;
    return result;
}

int main() {
    // Come Back in One Piece (★5)
    int N, M;
    cin >> N >> M;
    vector<vector<int>> graph(N);
    for (int i = 0; i < M; i++) {
        int a, b;
        cin >> a >> b;
        graph[a - 1].push_back(b - 1);
    }
    
    // メモ化用配列を初期化
    vector<vector<int>> memo(N, vector<int>(N, -1));

    // N頂点M辺の有向グラフ
    // 頂点(x,y)を行き来できる数を数える。
    int cnt = 0;
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            if (canreach(i, j, graph) && canreach(j, i, graph)) cnt++;
        }
    }
    cout << cnt << "\n";
}
#include <bits/stdc++.h>
using namespace std;

vector<vector<int>> graph, rgraph;
vector<bool> visited;
vector<int> order, comp;

void dfs1(int v) {
    visited[v] = true;
    for (int u : graph[v]) if (!visited[u]) dfs1(u);
    order.push_back(v);
}

void dfs2(int v, int c) {
    comp[v] = c;
    for (int u : rgraph[v]) if (comp[u] == -1) dfs2(u, c);
}

int main() {
    int N, M;
    cin >> N >> M;
    graph.resize(N);
    rgraph.resize(N);
    
    for (int i = 0; i < M; i++) {
        int a, b;
        cin >> a >> b;
        graph[a-1].push_back(b-1);
        rgraph[b-1].push_back(a-1);  // 逆辺
    }
    
    // Step 1: 帰りがけ順を記録
    visited.assign(N, false);
    for (int i = 0; i < N; i++) if (!visited[i]) dfs1(i);
    
    // Step 2: 逆グラフで帰りがけ順の逆順にDFS
    comp.assign(N, -1);
    int numSCC = 0;
    for (int i = N - 1; i >= 0; i--) {
        int v = order[i];
        if (comp[v] == -1) dfs2(v, numSCC++);
    }
    
    // Step 3: 各SCCのサイズを数えてペア数を計算
    vector<long long> sccSize(numSCC, 0);
    for (int i = 0; i < N; i++) sccSize[comp[i]]++;
    
    long long ans = 0;
    for (int i = 0; i < numSCC; i++) {
        ans += sccSize[i] * (sccSize[i] - 1) / 2;
    }
    cout << ans << "\n";
}
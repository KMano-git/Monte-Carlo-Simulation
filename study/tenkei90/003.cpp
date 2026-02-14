#include <bits/stdc++.h>
using namespace std;

// DFSで最も遠い頂点とその距離を返す
pair<int, int> dfs(int pos, int parent, int depth, vector<vector<int>> &adj) {
    int farthest_node = pos;
    int max_depth = depth;
    
    for (int next : adj[pos]) { // adj[pos]はposから行ける頂点のリスト
        if (next == parent) continue; // 来た道を戻らない
        
        auto [node, d] = dfs(next, pos, depth + 1, adj);
        if (d > max_depth) {
            max_depth = d;
            farthest_node = node;
        }
    }
    
    return {farthest_node, max_depth};
}

int main() {
    // Longest Circular Road (★4)
    // 木の直径を求める問題
    int N;
    cin >> N;
    
    // グラフの隣接リストを作成
    vector<vector<int>> adj(N);
    for (int i = 0; i < N - 1; i++) {  // 木なのでN-1本の辺
        int a, b;
        cin >> a >> b;
        a--; b--;  // 0-indexed に変換
        adj[a].push_back(b);
        adj[b].push_back(a);
    }

    // ステップ1: 頂点0から最も遠い頂点を見つける
    auto [node1, dist1] = dfs(0, -1, 0, adj);
    
    // ステップ2: その頂点から最も遠い頂点を見つける（これが直径）
    auto [node2, diameter] = dfs(node1, -1, 0, adj);
    
    // 答えは直径 + 1（新しい道を1本追加できるので）
    cout << diameter + 1 << endl;
    
    return 0;
}
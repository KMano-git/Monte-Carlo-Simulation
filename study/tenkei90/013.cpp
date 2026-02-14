#include <bits/stdc++.h>
using namespace std;

int main() {
    int N, M;
    cin >> N >> M;
    
    // 重み付き無向グラフ
    vector<vector<pair<int, int>>> graph(N);
    for (int i = 0; i < M; i++) {
        int a, b, c;
        cin >> a >> b >> c;
        a--;
        b--;
        graph.at(a).push_back({b, c});
        graph.at(b).push_back({a, c});
    }

    // ダイクストラ法
    auto dijkstra = [&](int start) -> vector<long long> {
        const long long INF = 1e18;
        vector<long long> dist(N, INF);
        // {距離, ノード} の優先度付きキュー（小さい順）
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> pq;
        
        dist[start] = 0;
        pq.push({0, start});
        
        while (!pq.empty()) {
            auto [d, v] = pq.top();
            pq.pop();
            
            // 既により短い経路が見つかっていればスキップ
            if (d > dist[v]) continue;
            
            // 隣接ノードを探索
            for (auto [next, cost] : graph[v]) {
                if (dist[v] + cost < dist[next]) {
                    dist[next] = dist[v] + cost;
                    pq.push({dist[next], next});
                }
            }
        }
        return dist;
    };
    
    // 街0からの最短距離
    vector<long long> dist_from_0 = dijkstra(0);
    // 街N-1からの最短距離（無向グラフなので、街iから街N-1への距離と同じ）
    vector<long long> dist_from_N1 = dijkstra(N - 1);
    
    // 各街iを経由する場合の最短経路長を出力
    for (int i = 0; i < N; i++) {
        cout << dist_from_0[i] + dist_from_N1[i] << "\n";
    }
}
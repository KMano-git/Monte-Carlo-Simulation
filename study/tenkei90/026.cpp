#include <bits/stdc++.h>
using namespace std;

int main() {
    int N;
    cin >> N;
    vector<vector<int>> graph(N);
    // 辺の本数はN-1本
    for (int i = 0; i < N - 1; i++) {
        int a,b;
        cin >> a >> b;
        a--;
        b--;
        graph.at(a).push_back(b);
        graph.at(b).push_back(a);
    }

    // N=2kのとき、隣り合わないk個の頂点を取り出す。
    // 木は二部グラフなので、二色に塗り分けることができる。
    // どちらかの色の頂点数は必ずN/2以上になる。
    // その色の頂点集合からk個選べばよい。

    vector<int> col(N, -1);
    vector<int> c0, c1;
    
    // BFSで彩色
    queue<int> q;
    q.push(0); // 始点
    col.at(0) = 0;
    c0.push_back(0); // 1-indexedなら+1するが、ここでは0-indexedで扱い、出力時に+1する

    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (int v : graph.at(u)) {
            if (col.at(v) == -1) {
                col.at(v) = 1 - col.at(u);
                if (col.at(v) == 0) c0.push_back(v);
                else c1.push_back(v);
                q.push(v);
            }
        }
    }

    vector<int> ans = c0;
    if (c1.size() > c0.size()) {
        ans = c1;
    }

    for (int i = 0; i < N / 2; i++) {
        cout << ans.at(i) + 1 << (i == N / 2 - 1 ? "" : " ");
    }
    cout << endl;
}
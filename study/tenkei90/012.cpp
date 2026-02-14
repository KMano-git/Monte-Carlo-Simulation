#include <bits/stdc++.h>
using namespace std;

class UnionFind {
public:
    vector<int> parent, rank_;
    
    UnionFind(int n) : parent(n), rank_(n, 0) {
        for (int i = 0; i < n; i++) parent[i] = i;
    }
    
    int find(int x) {
        if (parent[x] != x) parent[x] = find(parent[x]);
        return parent[x];
    }
    
    void unite(int x, int y) {
        int px = find(x), py = find(y);
        if (px == py) return;
        if (rank_[px] < rank_[py]) swap(px, py);
        parent[py] = px;
        if (rank_[px] == rank_[py]) rank_[px]++;
    }
    
    bool same(int x, int y) { return find(x) == find(y); }
};

int main() {
    int H, W;
    cin >> H >> W;
    
    UnionFind uf(H * W);
    vector<vector<bool>> grid(H, vector<bool>(W, false));
    
    int Q;
    cin >> Q;
    
    int dr[] = {-1, 1, 0, 0};
    int dc[] = {0, 0, -1, 1};
    
    while (Q--) {
        int t;
        cin >> t;
        if (t == 1) {
            int r, c;
            cin >> r >> c;
            r--; c--;
            grid[r][c] = true;
            // 隣接する赤いセルと結合
            for (int i = 0; i < 4; i++) {
                int nr = r + dr[i], nc = c + dc[i];
                if (nr >= 0 && nr < H && nc >= 0 && nc < W && grid[nr][nc]) {
                    uf.unite(r * W + c, nr * W + nc);
                }
            }
        } else {
            int r1, c1, r2, c2;
            cin >> r1 >> c1 >> r2 >> c2;
            r1--; c1--; r2--; c2--;
            if (grid[r1][c1] && grid[r2][c2] && uf.same(r1 * W + c1, r2 * W + c2)) {
                cout << "Yes\n";
            } else {
                cout << "No\n";
            }
        }
    }
}
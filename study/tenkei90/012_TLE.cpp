#include <bits/stdc++.h>
using namespace std;

bool canreach(int r1, int c1, int r2, int c2, vector<vector<bool>> &grid, vector<vector<bool>> &visited) {
    // 目標に到達したら true
    if (r1 == r2 && c1 == c2) return true;
    
    // 現在のセルを訪問済みにマーク
    visited.at(r1).at(c1) = true;
    
    // 4方向を探索（訪問済みでない赤いセルのみ）
    int dr[] = {-1, 1, 0, 0};
    int dc[] = {0, 0, -1, 1};
    
    for (int i = 0; i < 4; i++) {
        int nr = r1 + dr[i];
        int nc = c1 + dc[i];
        // 範囲内かつ赤いセルかつ未訪問の場合のみ探索
        int h = grid.size();
        int w = grid.at(0).size();
        if (nr >= 0 && nr < h && nc >= 0 && nc < w) {
            if (grid.at(nr).at(nc) && !visited.at(nr).at(nc)) {
                if (canreach(nr, nc, r2, c2, grid, visited)) {
                    return true;
                }
            }
        }
    }
    return false;
}

int main() {
    // Red Painting
    int H, W;
    cin >> H >> W;
    vector<vector<bool>> grid(H,vector<bool>(W,false));
    
    int Q;
    cin >> Q;
    for (int i = 0; i < Q; i++) {
        int t;
        cin >> t;
        if (t == 1) {
            int r, c;
            cin >> r >> c;
            grid.at(r - 1).at(c - 1) = true;
        } else {
            int r1, c1, r2, c2;
            cin >> r1 >> c1 >> r2 >> c2;
            if (grid.at(r1-1).at(c1-1) && grid.at(r2-1).at(c2-1)) {
                // 探索ごとにvisited配列を初期化
                vector<vector<bool>> visited(H, vector<bool>(W, false));
                if (canreach(r1-1, c1-1, r2-1, c2-1, grid, visited)) {
                    cout << "Yes" << endl;
                } else {
                    cout << "No" << endl;
                }
            } else {
                cout << "No" << endl;
            }
        }
    }

}
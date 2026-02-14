#include <bits/stdc++.h>
using namespace std;

int main() {
    int N;
    cin >> N;
    // 2次元いもす法 (2D Cumulative Sum)
    vector<vector<int>> grid(1010, vector<int>(1010, 0));
    for (int i = 0; i < N; i++) {
        int lx, ly, rx, ry;
        cin >> lx >> ly >> rx >> ry;
        // 加算と減算
        grid[lx][ly]++;
        grid[rx][ry]++;
        grid[lx][ry]--;
        grid[rx][ly]--;
    }

    // 横方向への累積和
    for (int i = 0; i <= 1000; i++) {
        for (int j = 1; j <= 1000; j++) {
            grid[i][j] += grid[i][j - 1];
        }
    }

    // 縦方向への累積和
    for (int j = 0; j <= 1000; j++) {
        for (int i = 1; i <= 1000; i++) {
            grid[i][j] += grid[i - 1][j];
        }
    }

    vector<int> ans(N + 1, 0);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            int k = grid[i][j];
            if (k > 0) ans[k]++;
        }
    }

    for (int k = 1; k <= N; k++) {
        cout << ans[k] << endl;
    }
}
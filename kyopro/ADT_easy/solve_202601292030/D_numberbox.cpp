#include <bits/stdc++.h>
using namespace std;

long long search(int n, int x, int y, int dx, int dy, vector<vector<int>> &A) {
    long long res = 0;
    int cur_x = x;
    int cur_y = y;
    for (int i = 0; i < n; i++) {
        if (cur_x == -1) {
            cur_x = n - 1;
        }
        if (cur_y == -1) {
            cur_y = n - 1;
        }
        if (cur_x == n) {
            cur_x = 0;
        }
        if (cur_y == n) {
            cur_y = 0;
        }
        res = res * 10 + A.at(cur_x).at(cur_y);
        cur_x += dx;
        cur_y += dy;
    }
    return res;
}

int main() {
    int n;
    cin >> n;
    vector<vector<int>> A(n, vector<int>(n));
    for(int i = 0; i < n; i++) {
        string s;
        cin >> s;
        for(int j = 0; j < n; j++) {
            A.at(i).at(j) = s[j] - '0';  // 文字を数字に変換
        }
    }

    // 一番大きい数字から向かう方向を決定してひたすら進む。
    // 8方向分はループする
    // 移動を定義、di,dj(x,y座標と順番が逆担っているので注意)
    pair<int, int> d[8] = {{1, 1}, {1, 0}, {1, -1}, {0, 1}, {0, -1}, {-1, 1}, {-1, 0}, {-1, -1}};

    long long max_num = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < 8; k++) {
                max_num = max(max_num, search(n, i, j, d[k].first, d[k].second, A));
            }
        }
    }
    cout << max_num << endl;
    

}
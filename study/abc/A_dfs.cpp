#include <bits/stdc++.h>
using namespace std;

void dfs(int x, int y, int h, int w, vector<vector<char>> &A, vector<vector<bool>> &reached) {
    //ゴール判定
    if (A.at(y).at(x) == 'g') {
        cout << "Yes" << endl;
        exit(0);
    }
    
    //壁判定
    if (A.at(y).at(x) == '#') {
        return;
    }
    
    //訪問済み判定
    if (reached.at(y).at(x)) {
        return;
    }
    
    //訪問済みマーク
    reached.at(y).at(x) = true;
    
    //上下左右の探索
    if (x + 1 < w) {
        dfs(x + 1, y, h, w, A, reached);
    }
    if (x - 1 >= 0) {
        dfs(x - 1, y, h, w, A, reached);
    }
    if (y + 1 < h) {
        dfs(x, y + 1, h, w, A, reached);
    }
    if (y - 1 >= 0) {
        dfs(x, y - 1, h, w, A, reached);
    }
}

int main() {
    //入力
    int h, w;
    cin >> h >> w;
    vector<vector<char>> A(h, vector<char>(w));
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            cin >> A.at(i).at(j);
        }
    }
    
    //dfs:深さ優先探索
    //startを探す
    int start_x = 0, start_y = 0;
    bool start_found = false;
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            if (A.at(i).at(j) == 's') {
                start_x = j;
                start_y = i;
                start_found = true;
                break;
            }
        }
        if (start_found) {
            break;
        }
    }
    
    //dfs
    vector<vector<bool>> reached(h, vector<bool>(w, false));
    dfs(start_x, start_y, h, w, A, reached);

    //ゴールに到達できなかった場合
    int goal_x = 0, goal_y = 0;
    bool goal_found = false;
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            if (A.at(i).at(j) == 'g') {
                goal_x = j;
                goal_y = i;
                goal_found = true;
                break;
            }
        }
        if (goal_found) {
            break;
        }
    }
    if (!reached.at(goal_y).at(goal_x)) {
        cout << "No" << endl;
    }
    
}
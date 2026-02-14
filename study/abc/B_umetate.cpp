#include <bits/stdc++.h>
using namespace std;

// A_dfs.cppを参考に
void dfs(int x, int y, int h, int w, vector<vector<char>> &A, vector<vector<bool>> &reached) {
    //海判定
    if (A.at(y).at(x) == 'x') {
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

bool is_ok(vector<vector<char>> &A) {
    
    // 陸地マスを1つ探す
    int y = 0, x = 0;
    bool found = false;
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            if (A.at(i).at(j) == 'o') {
                y = i;
                x = j;
                found = true;
                break;
            }
        }
        if (found) {
            break;
        }
    }
    
    //dfs
    vector<vector<bool>> reached(10, vector<bool>(10, false));
    dfs(x, y, 10, 10, A, reached);

    //陸地マスy,xから全ての陸地マスへ到達可能か判定
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            if (A.at(i).at(j) == 'o' && !reached.at(i).at(j)) {
                return false;
            }
        }
    }
    return true;
}

int main() {
    // 10x10のマス目の入力を受け取る
    vector<vector<char>> A(10, vector<char>(10));
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            cin >> A.at(i).at(j);
        }
    }

    bool ans = false;
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            if (A.at(i).at(j) == 'x') {
                A.at(i).at(j) = 'o';  // 埋め立て
                if (is_ok(A)) {
                    ans = true;
                }
                A.at(i).at(j) = 'x';  // 元に戻す
            }
        }
    }
    if (ans) {
        cout << "YES" << endl;
    } else {
        cout << "NO" << endl;
    }
}
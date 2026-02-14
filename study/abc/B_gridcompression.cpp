#include <bits/stdc++.h>
using namespace std;

int main() {
    // input
    int h, w;
    cin >> h >> w;
    vector<vector<char>> a(h, vector<char>(w));
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            cin >> a.at(i).at(j);
        }
    }

    // compression
    vector<bool> iswhite_row(h, true);
    vector<bool> iswhite_col(w, true);
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            if (a.at(i).at(j) == '#') {
                iswhite_row.at(i) = false;
                iswhite_col.at(j) = false;
            }
        }
    }
    // 行の削除（後ろから削除すればインデックスがずれない）
    for (int i = h - 1; i >= 0; i--) {
        if (iswhite_row.at(i)) {
            a.erase(a.begin() + i);
        }
    }
    
    // 列の削除（後ろから削除）
    for (int j = w - 1; j >= 0; j--) {
        if (iswhite_col.at(j)) {
            for (int i = 0; i < (int)a.size(); i++) {
                a.at(i).erase(a.at(i).begin() + j);
            }
        }
    }

    // output（a.size() を使う）
    for (int i = 0; i < (int)a.size(); i++) {
        for (int j = 0; j < (int)a.at(i).size(); j++) {
            cout << a.at(i).at(j);
        }
        cout << endl;
    }
}
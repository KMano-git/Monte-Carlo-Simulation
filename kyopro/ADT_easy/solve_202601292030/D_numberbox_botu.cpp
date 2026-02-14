#include <bits/stdc++.h>
using namespace std;

int max_8dir(int a1, int a2, int a3, int a4, int a5, int a6, int a7, int a8) {
    int max_val = a1;
    if (a2 > max_val) max_val = a2;
    if (a3 > max_val) max_val = a3;
    if (a4 > max_val) max_val = a4;
    if (a5 > max_val) max_val = a5;
    if (a6 > max_val) max_val = a6;
    if (a7 > max_val) max_val = a7;
    if (a8 > max_val) max_val = a8;
    return max_val;
}

int main() {
    int n;
    cin >> n;
    vector<vector<int>> A(n, vector<int>(n));
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            cin >> A.at(i).at(j);
        }
    }

    // ループつきの探索 // 大きいマスを選べばいいので探索ではないか
    // 一番大きい数字のマスを見つける
    int max_val = 0;
    int max_x = 0;
    int max_y = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (A.at(i).at(j) > max_val) {
                max_val = A.at(i).at(j);
                max_x = i;
                max_y = j;
            }
        }
    }

    // 一番大きい数字の候補を列挙する
    vector<pair<int, int>> candidates;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (A.at(i).at(j) == max_val) {
                candidates.push_back({i, j});
            }
        }
    }

    // 候補のマスからスタートしてnマス分探索する
    // 8方向分ループしてもいいようにマスを拡張する
    vector<vector<int>> A_ext(3 * n, vector<int>(3 * n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A_ext.at(i + n).at(j + n) = A.at(i).at(j);
        }
    }
    
    int number_of_candidates = candidates.size();
    vector<int> max_scores(number_of_candidates, 0);
    for (int i = 0; i < number_of_candidates; i++) {
        int x = candidates.at(i).first + n;
        int y = candidates.at(i).second + n;
        int current_score = A_ext.at(x).at(y) * pow(10, n - 1);
        for (int j = 0; j < n; j++) {
            current_score += max_8dir(A_ext.at(x - 1).at(y - 1), A_ext.at(x - 1).at(y), A_ext.at(x - 1).at(y + 1), A_ext.at(x).at(y - 1), A_ext.at(x).at(y + 1), A_ext.at(x + 1).at(y - 1), A_ext.at(x + 1).at(y), A_ext.at(x + 1).at(y + 1)) * pow(10, n - 2 - j);
            // これじゃあ、探索が続けられないのでボツ
        }
        max_scores.at(i) = current_score;
    }

    sort(max_scores.rbegin(), max_scores.rend());
    cout << max_scores.at(0) << endl;
    
}
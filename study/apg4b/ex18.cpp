#include <bits/stdc++.h>
using namespace std;

int main() {
    // 入力受付
    int n, m;
    cin >> n >> m;
    vector<int> a(m), b(m);
    for (int i = 0; i < m; i++) {
        cin >> a.at(i) >> b.at(i);
    }
    // 試合結果の表を作る
    vector<vector<char>> result(n, vector<char>(n,'-'));
    for (int i = 0; i < m; i++) {
        result.at(a.at(i)-1).at(b.at(i)-1) = 'o';
        result.at(b.at(i)-1).at(a.at(i)-1) = 'x';
    }
    // 結果を出力
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j == n-1) {
                cout << result.at(i).at(j) << endl;
            } else {
                cout << result.at(i).at(j) << " ";
            }
        }
    }
}
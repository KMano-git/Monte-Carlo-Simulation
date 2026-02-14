#include <bits/stdc++.h>
using namespace std;

int main() {
    int h, w;
    cin >> h >> w;
    // 2次元配列を入力から受け取る
    vector<vector<char>> s(h, vector<char>(w));
    for (int i = 0; i < h; i++){
        for (int j = 0; j < w; j++){
            cin >> s[i][j];
        }
    }
    vector<vector<char>> ans(h, vector<char>(w));
    for (int i = 0; i < h; i++){
        for (int j = 0; j < w; j++){
            int count = 0;
            if (s[i][j] == '#'){
                ans.at(i).at(j) = '#';
            } else {
                for (int di = -1; di <= 1; di++){
                    for (int dj = -1; dj <= 1; dj++){
                        if (di == 0 && dj == 0){
                            continue;
                        }
                        int ni = i + di;
                        int nj = j + dj;
                        if (ni < 0 || ni >= h || nj < 0 || nj >= w){
                            continue;
                        }
                        if (s.at(ni).at(nj) == '#'){
                            count++;
                        }
                    }
                }
                ans.at(i).at(j) = count + '0';
            }
        }
    }
    // 出力
    for (int i = 0; i < h; i++){
        for (int j = 0; j < w; j++){
            cout << ans.at(i).at(j);
        }
        cout << endl;
    }
}
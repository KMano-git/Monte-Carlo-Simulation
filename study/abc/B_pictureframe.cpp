#include <bits/stdc++.h>
using namespace std;

int main() {
    int h, w;
    cin >> h >> w;
    vector<vector<char>> a(h, vector<char>(w));
    for (int i = 0; i < h; i++){
        for (int j = 0; j < w; j++){
            cin >> a.at(i).at(j);
        }
    }
    vector<vector<char>> b(h+2, vector<char>(w+2));
    for (int i = 0; i < h+2; i++){
        for (int j = 0; j < w+2; j++) {
            if (i == 0 || i == h+1 || j == 0 || j == w+1){
                b.at(i).at(j) = '#';
            } else {
                b.at(i).at(j) = a.at(i-1).at(j-1);
            }
        }
    }
    for (int i = 0; i < h+2; i++){
        for (int j = 0; j < w+2; j++){
            cout << b.at(i).at(j);
        }
        cout << endl;
    }

}
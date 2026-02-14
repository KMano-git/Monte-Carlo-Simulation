#include <bits/stdc++.h>
using namespace std;

int main() {
    // Cross Sum (★2)
    int H, W;
    cin >> H >> W;
    vector<vector<int>> A(H, vector<int>(W));
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            cin >> A.at(i).at(j);
        }
    }
    
    vector<int> row_sum(H, 0);
    vector<int> col_sum(W, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            row_sum.at(i) += A.at(i).at(j);
            col_sum.at(j) += A.at(i).at(j);
        }
    }
    
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            cout << row_sum.at(i) + col_sum.at(j) - A.at(i).at(j) << " ";
        }
        cout << endl;
    }
}
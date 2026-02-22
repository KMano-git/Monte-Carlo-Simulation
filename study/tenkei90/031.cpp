#include <bits/stdc++.h>
using namespace std;

int G[55][1550];
int mex[1550];

void init() {
    for (int w = 0; w <= 50; w++) {
        for (int b = 0; b <= 1500; b++) {
            if (w >= 1 && b + w <= 1500) {
                mex[G[w - 1][b + w]] = 1;
            }
            if (b >= 2) {
                for (int k = 1; k <= b / 2; k++) {
                    mex[G[w][b - k]] = 1;
                }
            }
            int g = 0;
            while (mex[g]) g++;
            G[w][b] = g;
            
            if (w >= 1 && b + w <= 1500) {
                mex[G[w - 1][b + w]] = 0;
            }
            if (b >= 2) {
                for (int k = 1; k <= b / 2; k++) {
                    mex[G[w][b - k]] = 0;
                }
            }
        }
    }
}

int main() {
    int N;
    cin >> N;

    vector<int> W(N), B(N);
    for (int i = 0; i < N; i++) {
        cin >> W.at(i);
    }
    for (int i = 0; i < N; i++) {
        cin >> B.at(i);
    }

    init();

    int xor_sum = 0;
    for (int i = 0; i < N; i++) {
        xor_sum ^= G[W[i]][B[i]];
    }

    if (xor_sum != 0) {
        cout << "First" << endl;
    } else {
        cout << "Second" << endl;
    }

    return 0;
}
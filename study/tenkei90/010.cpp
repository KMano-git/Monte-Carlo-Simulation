#include <bits/stdc++.h>
using namespace std;

int main() {
    int N;
    cin >> N;
    vector<int> C(N), P(N); // Ci = 1 or 2, Pi = score
    for (int i = 0; i < N; i++) {
        cin >> C[i] >> P[i];
    }

    int Q;
    cin >> Q;
    vector<int> L(Q), R(Q);
    for (int i = 0; i < Q; i++) {
        cin >> L[i] >> R[i];
    }

    // 計算量削減のために累積和を使う
    vector<int> sum1(N + 1, 0);
    vector<int> sum2(N + 1, 0);
    for (int i = 0; i < N; i++) {
        if (C.at(i) == 1) {
            sum1.at(i + 1) = sum1.at(i) + P.at(i);
            sum2.at(i + 1) = sum2.at(i);
        } else {
            sum1.at(i + 1) = sum1.at(i);
            sum2.at(i + 1) = sum2.at(i) + P.at(i);
        }
    }
    
    for (int i = 0; i < Q; i++) {
        int score1 = sum1.at(R.at(i)) - sum1.at(L.at(i) - 1);
        int score2 = sum2.at(R.at(i)) - sum2.at(L.at(i) - 1);
        cout << score1 << " " << score2 << endl;
    }
}
#include <bits/stdc++.h>
using namespace std;

int main() {
    int N, K;
    cin >> N >> K;
    vector<int> A(N), B(N);
    for(int i = 0; i < N; i++) cin >> A[i];
    for(int i = 0; i < N; i++) cin >> B[i];
    
    int cnt = 0;
    for (int i = 0; i < N; i++) {
        cnt += abs(A.at(i) - B.at(i));
    }

    if (cnt > K) cout << "No" << "\n";
    else if ((K - cnt) % 2 == 1) cout << "No" << "\n";
    else cout << "Yes" << "\n";
}
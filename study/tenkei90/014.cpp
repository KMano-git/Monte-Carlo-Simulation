#include <bits/stdc++.h>
using namespace std;

int main() {
    int N;
    cin >> N;
    vector<int> A(N), B(N);
    for (int i = 0; i < N; i++) {
        cin >> A[i];
    }
    for (int i = 0; i < N; i++) {
        cin >> B[i];
    }
    
    // 組み合わせ＋最小値問題　→　貪欲法
    // Aをソート
    sort(A.begin(), A.end());
    // Bをソート
    sort(B.begin(), B.end());
    
    // A[i]とB[i]の差の絶対値を足していく
    long long ans = 0;
    for (int i = 0; i < N; i++) {
        ans += abs(A.at(i) - B.at(i));
    }
    
    cout << ans << endl;
}
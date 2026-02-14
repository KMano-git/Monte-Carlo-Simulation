#include <bits/stdc++.h>
using namespace std;

int main() {
    // CP classes (★3)
    int N;
    cin >> N;
    vector<int> A(N);
    for (int i = 0; i < N; i++) {
        cin >> A[i];
    }

    sort(A.begin(), A.end());
    int Q;
    cin >> Q;
    vector<int> B(Q);
    for (int i = 0; i < Q; i++) {
        cin >> B[i];
        auto it = lower_bound(A.begin(), A.end(), B[i]);
        int idx = it - A.begin();
        if (idx == N) {
            cout << B.at(i) - A.at(N - 1) << endl;
        } else if (idx == 0) {
            cout << A.at(0) - B.at(i) << endl;
        } else {
            cout << min(B.at(i) - A.at(idx-1), A.at(idx) - B.at(i)) << endl;
        }
    }
}
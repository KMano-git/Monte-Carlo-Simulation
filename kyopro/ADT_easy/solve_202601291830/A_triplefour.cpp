#include <bits/stdc++.h>
using namespace std;

int main() {
    // 3 <= n <= 100
    // 1 <= A[i] <= 100
    int n;
    cin >> n;
    vector<int> A(n);
    for (int i = 0; i < n; i++) {
        cin >> A.at(i);
    }

    bool ans = false;
    for (int i = 0; i < n-2; i++) {
        if (A.at(i) == A.at(i+1) && A.at(i) == A.at(i+2)) {
            ans = true;
            break;
        }
    }
    if (ans) cout << "Yes" << endl;
    else cout << "No" << endl;
}
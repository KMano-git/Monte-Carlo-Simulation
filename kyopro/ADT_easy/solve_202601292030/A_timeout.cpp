#include <bits/stdc++.h>
using namespace std;

int main() {
    int n , s;
    cin >> n >> s;
    vector<int> T(n);
    for(int i = 0; i < n; i++) {
        cin >> T.at(i);
    }
    
    if (T.at(0) > s) {
        cout << "No" << endl;
        return 0;
    }

    bool wake = false;
    for (int i = 0; i < n-1; i++) {
        if (T.at(i+1) - T.at(i) > s) {
            wake = true;
            break;
        }
    }

    if (wake) {
        cout << "No" << endl;
    } else {
        cout << "Yes" << endl;
    }
}
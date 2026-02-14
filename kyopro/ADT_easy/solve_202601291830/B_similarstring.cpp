#include <bits/stdc++.h>
using namespace std;

bool issimilar(char a, char b) {
    if (a == b) return true;
    if (a == '1' && b == 'l') return true;
    if (a == 'l' && b == '1') return true;
    if (a == '0' && b == 'o') return true;
    if (a == 'o' && b == '0') return true;
    return false;
}

int main() {
    int n;
    cin >> n;
    string s, t;
    cin >> s >> t;

    bool ans = true;
    for (int i = 0; i < n; i++) {
        if (!issimilar(s.at(i), t.at(i))) {
            ans = false;
            break;
        }
    }
    if (ans) cout << "Yes" << endl;
    else cout << "No" << endl;
}
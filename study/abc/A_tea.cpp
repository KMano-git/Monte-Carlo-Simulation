#include <bits/stdc++.h>
using namespace std;

int main() {
    string s;
    cin >> s;
    int n = s.size();
    if (s.at(n-1) == 'T') {
        cout << "YES" << endl;
    } else {
        cout << "NO" << endl;
    }
}
#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin >> n;
    int count = 0;

    for (int i = 0; i < n; i++) {
        string s;
        cin >> s;
        if (s == "Y") {
            count += 1;
        } else {
            continue;
        }
    }

    if (count >= 1) {
        cout << "Four" << endl;
    } else {
        cout << "Three" << endl;
    }
}
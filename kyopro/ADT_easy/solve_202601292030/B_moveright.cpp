#include <bits/stdc++.h>
using namespace std;

int main() {
    string s;
    cin >> s;

    for (int i = 3; i >= 1; i--) {
        s.at(i) = s.at(i - 1);
    }
    s.at(0) = '0';
    cout << s << endl;
}
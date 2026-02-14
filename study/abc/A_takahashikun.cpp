#include <bits/stdc++.h>
using namespace std;

int main() {
    string a, b;
    cin >> a >> b;
    int na = a.size();
    int nb = b.size();
    if (na > nb) {
        cout << a << endl;
    } else {
        cout << b << endl;
    }
}
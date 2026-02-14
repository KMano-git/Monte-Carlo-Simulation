#include <bits/stdc++.h>
using namespace std;

int main() {
    int a, n;
    cin >> n >> a;
    int amari = n % 500;
    if (amari <= a) {
        cout << "Yes" << endl;
    } else {
        cout << "No" << endl;
    }
}
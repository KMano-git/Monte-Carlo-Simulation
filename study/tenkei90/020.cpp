#include <bits/stdc++.h>
using namespace std;
using ll = long long;

ll power(ll a, ll b) {
    ll res = 1;
    for (int i = 0; i < b; i++) {
        res *= a;
    }
    return res;
}

int main() {
    ll a, b, c;
    cin >> a >> b >> c;

    if (power(c, b) > a) cout << "Yes" << "\n";
    else cout << "No" << "\n";
}
#include <bits/stdc++.h>
using namespace std;
using ll = long long;

ll power(ll x, ll y) {
    ll res = 1;
    for (int i = 0; i < y; i++) {
        res *= x;
    }
    return res;
}

int main() {
    ll k;
    cin >> k;
    ll a, b;
    cin >> a >> b;

    ll num_a = 0, num_b = 0;
    ll base = 1;

    // aをk進数として解釈
    while (a > 0) {
        num_a += (a % 10) * base;
        a /= 10;
        base *= k;
    }
    base = 1;
    // bをk進数として解釈
    while (b > 0) {
        num_b += (b % 10) * base;
        b /= 10;
        base *= k;
    }

    cout << num_a * num_b << endl;
}
// 問題解釈にミスがあったみたい。値として1<=a,b<=10^5なので、6桁以上になることもある。

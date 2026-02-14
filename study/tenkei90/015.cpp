#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const ll MOD = 1e9 + 7;

ll power(ll x, ll y) {
    ll res = 1;
    x %= MOD;
    while (y > 0) {
        if (y & 1) res = (res * x) % MOD;
        y >>= 1;
        x = (x * x) % MOD;
    }
    return res;
}

int main() {
    int N;
    cin >> N;
    
    // 階乗と逆元を前計算
    vector<ll> fact(N + 1), inv_fact(N + 1);
    fact[0] = 1;
    for (int i = 1; i <= N; i++) {
        fact[i] = fact[i - 1] * i % MOD;
    }
    inv_fact[N] = power(fact[N], MOD - 2);
    for (int i = N - 1; i >= 0; i--) {
        inv_fact[i] = inv_fact[i + 1] * (i + 1) % MOD;
    }
    
    // nCr を O(1) で計算
    auto nCr = [&](int n, int r) -> ll {
        if (r < 0 || r > n) return 0;
        return fact[n] * inv_fact[r] % MOD * inv_fact[n - r] % MOD;
    };
    
    // 各 k について答えを計算
    for (int k = 1; k <= N; k++) {
        ll ans = 0;
        for (int m = 1; N - (ll)(m - 1) * (k - 1) >= m; m++) {
            // N - (m-1)*(k-1) 個から m 個選ぶ
            ans = (ans + nCr(N - (m - 1) * (k - 1), m)) % MOD;
        }
        cout << ans << "\n";
    }
    
    return 0;
}
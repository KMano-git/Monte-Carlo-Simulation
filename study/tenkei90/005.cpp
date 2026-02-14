#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const ll MOD = 1e9 + 7;

map<pair<ll, ll>, ll> memo; // メモ化用のテーブル

// current_digit: 現在の桁数
// remainder: 現在の数をBで割った余り
ll solve(ll N, ll B, ll K, vector<int> &c, ll current_digit, ll remainder) {
    // ベースケース：N桁すべて決定した
    if (current_digit == N) {
        return (remainder == 0) ? 1 : 0;
    }
    
    // メモ化チェック
    auto key = make_pair(current_digit, remainder);
    if (memo.count(key)) {
        return memo[key];
    }
    
    ll cnt = 0;
    for (int i = 0; i < K; i++) {
        ll new_remainder = (remainder * 10 + c[i]) % B;
        cnt = (cnt + solve(N, B, K, c, current_digit + 1, new_remainder)) % MOD;
    }
    
    // 結果をメモに保存
    memo[key] = cnt;
    return cnt;
}

int main() {
    ll N, B, K;
    cin >> N >> B >> K;
    vector<int> c(K);
    for (int i = 0; i < K; i++) {
        cin >> c[i];
    }
    
    cout << solve(N, B, K, c, 0, 0) << endl;
    return 0;
}
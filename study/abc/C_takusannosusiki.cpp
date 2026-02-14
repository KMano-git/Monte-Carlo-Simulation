#include <bits/stdc++.h>
using namespace std;
using ll = long long;

ll ll_pow(ll base, ll exp) {
    ll result = 1;
    for (ll i = 0; i < exp; i++) {
        result *= base;
    }
    return result;
}

ll sum_allpattern(ll n, vector<ll> &a) {
    if (n == 1) { // ベースケース
        return a.at(0);
    }
    ll ans = 0;
    ans += 2 * sum_allpattern(n-1, a);
    // i桁分の先頭に付け足される分を考える
    for (ll i = 0; i < n-1; i++) {
        ans += ll_pow(10,i) * ll_pow(2, n-i-2) * a.at(n-1);
    }
    ans += ll_pow(10, n-1) * a.at(n-1);
    return ans;
}

int main() {
    // 入力
    string s;
    cin >> s;
    // 文字列の数字をllに変換
    ll n = s.size();
    for (ll i = 0; i < n; i++) {
        s.at(i) = s.at(i) - '0';
    }
    // 数字列を数字に変換
    ll num = 0;
    for (ll i = 0; i < n; i++) {
        num = num * 10 + s.at(i);
    }
    vector<ll> a(10); // a(i) = 10^i の位の数字
    for (ll i = 0; i < n; i++) {
        a.at(i) = s.at(n - i - 1);
    }
    
    // 再帰処理
    ll ans = 0;
    ans = sum_allpattern(n, a);

    // // 確かめ用
    // cout << num << endl;
    // for (ll i = 0; i < n; i++) {
    //     cout << a.at(i) << " ";
    // }
    // cout << endl;

    // 答えを出力
    cout << ans << endl;
}
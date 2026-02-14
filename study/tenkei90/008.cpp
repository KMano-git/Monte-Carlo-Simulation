#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const ll MOD = 1e9 + 7;

int main() {
    // AtCounter (★4)
    int N;
    cin >> N;
    string S;
    cin >> S;

    //文字列を抜き出して部分列"atcoder"の数を数える→動的計画法
    string T = "atcoder";
    vector<vector<ll>> dp(N + 1, vector<ll>(8, 0));
    
    for (int i = 0; i <= N; i++) {
        dp[i][0] = 1;
    }
    
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= 7; j++) {
            dp[i][j] = dp[i - 1][j];
            if (S[i - 1] == T[j - 1]) {
                dp[i][j] = (dp[i][j] + dp[i - 1][j - 1]) % MOD;
            }
        }
    }
    
    cout << dp[N][7] << endl;
}
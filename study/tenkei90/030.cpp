#include <bits/stdc++.h>
using namespace std;

int main() {
    int n, k;
    cin >> n >> k;
    
    // 2以上N以下の整数の内、K種類以上の素因数を持つものの個数を求める。
    // 2 <= N <= 1e7
    vector<int> c(n + 1, 0);
    int ans = 0;
    
    for (int i = 2; i <= n; ++i) {
        if (c.at(i) == 0) {
            for (int j = i; j <= n; j += i) {
                c.at(j)++;
            }
        }
        if (c.at(i) >= k) {
            ans++;
        }
    }
    
    cout << ans << endl;
    return 0;
}
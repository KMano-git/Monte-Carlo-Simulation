#include <bits/stdc++.h>
using namespace std;

int main() {
    int N;
    cin >> N;
    
    vector<tuple<int, int, int>> jobs(N);  // {deadline, cost, score}
    int D_max = 0;
    for (int i = 0; i < N; i++) {
        int d, c, s;
        cin >> d >> c >> s;
        jobs[i] = {d, c, s};
        D_max = max(D_max, d);
    }
    
    // 締め切り順にソート
    sort(jobs.begin(), jobs.end());
    
    // dp[j] = j日使ったときの最大報酬
    vector<long long> dp(D_max + 1, 0);
    
    for (int i = 0; i < N; i++) {
        auto [d, c, s] = jobs[i];
        
        // 逆順（0-1ナップサック）
        for (int j = d; j >= c; j--) {
            dp[j] = max(dp[j], dp[j - c] + s);
        }
    }
    
    cout << *max_element(dp.begin(), dp.end()) << endl;
}
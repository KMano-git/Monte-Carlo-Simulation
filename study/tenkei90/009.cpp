#include <bits/stdc++.h>
using namespace std;
using ll = long long;

int main() {
    ll N;
    cin >> N;
    vector<ll> x(N), y(N);
    for (ll i = 0; i < N; i++) {
        cin >> x[i] >> y[i];
    }

    double max_angle = 0;
    
    // 中心点 j を固定
    for (int j = 0; j < N; j++) {
        vector<double> theta;
        
        // j から各点への偏角を計算
        for (int i = 0; i < N; i++) {
            if (i == j) continue;
            theta.push_back(atan2(y[i] - y[j], x[i] - x[j]));
        }
        sort(theta.begin(), theta.end());
        
        // 周期性対応: 配列を2倍に
        int m = theta.size();
        for (int k = 0; k < m; k++) {
            theta.push_back(theta[k] + 2 * M_PI);
        }
        
        // 各角度に対して180°に最も近いペアを探す
        for (int k = 0; k < m; k++) {
            auto it = lower_bound(theta.begin(), theta.end(), theta[k] + M_PI);
            int idx = it - theta.begin();
            
            // idx と idx-1 の両方を確認
            double diff = theta[idx] - theta[k];  // π より大きいか等しい
            max_angle = max(max_angle, min(diff, 2 * M_PI - diff));
            
            if (idx > 0) {
                diff = theta[idx - 1] - theta[k];  // π より小さい
                max_angle = max(max_angle, min(diff, 2 * M_PI - diff));
            }
        }
    }
    
    // ラジアンから度に変換
    cout << fixed << setprecision(10) << max_angle * 180.0 / M_PI << endl;
}
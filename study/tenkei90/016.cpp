#include <bits/stdc++.h>
using namespace std;

int main() {
    long long N;
    cin >> N;
    long long A, B, C;
    cin >> A >> B >> C;

    // ちょうどN円をA,B,Cで払う
    long long ans = 1e9;
    
    for (int i = 0; i * A <= N && i < 1e4; i++) {
        for (int j = 0; j * B <= N - i * A && j < 1e4; j++) {
            long long rem = N - i * A - j * B;
            if (rem % C == 0 && rem >= 0 && rem / C < 1e4) {
                ans = min(ans, i + j + rem / C);
            }
        }
    }

    cout << ans << "\n";
}
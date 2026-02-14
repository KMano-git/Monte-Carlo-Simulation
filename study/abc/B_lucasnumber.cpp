#include <bits/stdc++.h>
using namespace std;

//再帰
// long long lucasnumber(long long n) {
//     if (n == 0) return 2;
//     if (n == 1) return 1;
//     return lucasnumber(n - 1) + lucasnumber(n - 2);
// }
// int main() {
//     long long n;
//     cin >> n;

//     cout << lucasnumber(n) << endl;
// }

// 一般項を使用

// int main() {
//     long long n;
//     cin >> n;

//     long long lucasnumber;
//     double alpha = (1 + sqrt(5)) / 2;
//     double beta = (1 - sqrt(5)) / 2;
//     lucasnumber = pow(beta, n) + pow(alpha, n);
//     cout << lucasnumber << endl;
// }

//ループを使用
int main() {
    long long n;
    cin >> n;

    if (n == 0) {
        cout << 2 << endl;
        return 0;
    }
    if (n == 1) {
        cout << 1 << endl;
        return 0;
    }

    // ループで計算（メモ化と同じ効率）
    long long prev2 = 2;  // L(0)
    long long prev1 = 1;  // L(1)
    long long current;
    
    for (long long i = 2; i <= n; i++) {
        current = prev1 + prev2;
        prev2 = prev1;
        prev1 = current;
    }
    
    cout << current << endl;
}
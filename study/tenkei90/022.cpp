#include <bits/stdc++.h>
using namespace std;
using ll = long long;

int main() {
    ll A, B, C;
    cin >> A >> B >> C;
    
    // A,B,Cの最大公約数を求める
    ll g = __gcd(A, B);
    g = __gcd(g, C);
    
    // gの大きさにA,B,Cをそれぞれ切る
    ll cut = A / g + B / g + C / g - 3;
    cout << cut << "\n";
}
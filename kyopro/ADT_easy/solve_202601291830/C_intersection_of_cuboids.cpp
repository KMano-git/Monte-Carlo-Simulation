#include <bits/stdc++.h>
using namespace std;

int main() {
    int a, b, c, d, e, f;
    cin >> a >> b >> c >> d >> e >> f;
    int g, h, i, j, k, l;
    cin >> g >> h >> i >> j >> k >> l;

    // 直方体(a,b,c,d,e,f)と(g,h,i,j,k,l)の重なっている部分の面積を考えたい
    // どちらもxy,yz,zx平面のいずれにも平行な直方体
    // (a<x<d, b<y<e, c<z<f) と (g<x<j, h<y<k, i<z<l) の重なっている部分の面積を考えたい
    
    // 重なっているか判定；
    bool iscollision = false;
    if ( (j - a) * (d - g) > 0 && (k - b) * (e - h) > 0 && (l - c) * (f - i) > 0) iscollision = true;
    if (iscollision) cout << "Yes" << endl;
    else cout << "No" << endl;
}
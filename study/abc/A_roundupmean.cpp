#include <bits/stdc++.h>
using namespace std;

int main() {
    int a,b;
    cin >> a >> b;
    double x = (a+b)/2.0;
    cout << ceil(x) << endl;
    // (a + b - 1) / b で a/b の切り上げ
    // int x = (a + b - 1) / b;  
    // 整数演算のみで切り上げ除算
}
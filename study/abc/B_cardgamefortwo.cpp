#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin >> n;
    vector<int> a(n);
    for (int i = 0; i < n; i++) {
        cin >> a.at(i);
    }
    
    // 降順にソート
    sort(a.begin(), a.end(), greater<int>());
    
    int a_point = 0;
    int b_point = 0;
    for (int i = 0; i < n; i++) {
        if (i % 2 == 0) {
            a_point += a.at(i);  // 偶数番目（0,2,4...）はAの手番
        } else {
            b_point += a.at(i);  // 奇数番目（1,3,5...）はBの手番
        }
    }
    cout << a_point - b_point << endl;
}
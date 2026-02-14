#include <bits/stdc++.h>
using namespace std;

int main() {
    int n, m;
    cin >> n >> m;
    vector<int> a(m), b(m);
    for (int i = 0; i < m; i++) {
        cin >> a.at(i) >> b.at(i);
    }
    //問題設定がよくわからなかったが、aiとbiの間に道があるということらしい
    //つまり、数字が出てきた数を数えるだけでいい
    vector<int> city(n,0);
    for (int i = 0; i < m; i++) {
        city.at(a.at(i)-1)++;
        city.at(b.at(i)-1)++;
    }
    for (int i = 0; i < n; i++) {
        cout << city.at(i) << endl;
    }
}
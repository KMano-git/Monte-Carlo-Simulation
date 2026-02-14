#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin >> n;
    vector<string> s(n);
    for (int i = 0; i < n; i++) {
        cin >> s.at(i);
    }
    int m;
    cin >> m;
    vector<string> t(m);
    for (int i = 0; i < m; i++) {
        cin >> t.at(i);
    }

    // sの内、どの文字列を選択すれば良いかを数える。
    int maxCount = 0;
    for (int i = 0; i < n; i++) {
        int count = 0;
        for (int j = 0; j < n; j++) {
            if (s.at(i) == s.at(j)) {
                count++;
            }
        }
        for (int k = 0; k < m; k++) {
            if (s.at(i) == t.at(k)) {
                count--;
            }
        }
        maxCount = max(maxCount, count);
    }
    cout << maxCount << endl;
}
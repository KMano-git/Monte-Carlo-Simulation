#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin >> n;
    vector<int> a(n);
    for (int i = 0; i < n; i++) {
        cin >> a.at(i);
    }
    
    int count = 0;
    bool flag = true;
    while (flag) {
        for (int i = 0; i < n; i++) {
            if (a.at(i) % 2 != 0) {
                flag = false;
                break;
            }
        }
        if (!flag) break;
        for (int i = 0; i < n; i++) {
            a.at(i) = a.at(i) / 2;
        }
        count++;
    }
    cout << count << endl;
}
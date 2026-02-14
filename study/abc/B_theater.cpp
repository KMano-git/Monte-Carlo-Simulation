#include <bits/stdc++.h>
using namespace std;

int main() {
    int n, l, r, count = 0;
    cin >> n;
    for (int i = 0; i < n; i++) {
        cin >> l >> r;
        count += r - l + 1;
    }
    cout << count << endl;
}
#include <bits/stdc++.h>
using namespace std;

int main() {
    int a, b, c;
    cin >> a >> b >> c;
    int x;
    cin >> x; // xは50の倍数

    int count = 0;
    int target = x / 50;
    for (int i = 0; i <= a; i++) {
        for (int j = 0; j <= b; j++) {
            for (int k = 0; k <= c; k++) {
                if (k + 2 * j + 10 * i == target) count++;
            }
        }
    }
    cout << count << endl;
}
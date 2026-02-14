#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin >> n;
    int count_min = 1000;
    for (int i = 0; i < n; i++) {
        int a;
        cin >> a;
        int count = 0;
        while (a % 2 == 0) {
            a /= 2;
            count++;
        }
        if (count < count_min) {
            count_min = count;
        }
    }
    cout << count_min << endl;
}
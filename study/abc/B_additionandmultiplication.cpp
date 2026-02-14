#include <bits/stdc++.h>
using namespace std;

int main() {
    int n, k;
    cin >> n >> k;
    int number = 1;
    for (int i = 0; i < n; i++) {
        if (number >= k) {
            number += k;   
        } else {
            number *= 2;
        }
    }
    cout << number << endl;
}
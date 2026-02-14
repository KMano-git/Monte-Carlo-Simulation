#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin >> n;
    int sum = 0;
    sum += n % 10;
    sum += n / 10 % 10;
    sum += n / 100 % 10;
    sum += n / 1000 % 10;
    sum += n / 10000 % 10;
    sum += n / 100000 % 10;
    sum += n / 1000000 % 10;
    sum += n / 10000000 % 10;
    sum += n / 100000000 % 10;
    if (n % sum == 0) {
        cout << "Yes" << endl;
    } else {
        cout << "No" << endl;
    }
}
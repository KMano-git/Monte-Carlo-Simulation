#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin >> n;
    int s1, s2, s3, s4;
    s1 = n / 1000;
    s2 = (n / 100) % 10;
    s3 = (n / 10) % 10;
    s4 = n % 10;
    if (s1 == s2 && s2 == s3) {
        cout << "Yes" << endl;
    } else if (s2 == s3 && s3 == s4) {
        cout << "Yes" << endl;
    } else {
        cout << "No" << endl;
    }
}
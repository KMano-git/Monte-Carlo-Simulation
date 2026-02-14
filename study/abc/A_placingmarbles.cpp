#include <bits/stdc++.h>
using namespace std;

int main() {
    int s1, s2, s3;
    int input_s;
    cin >> input_s;
    s1 = input_s / 100;
    s2 = (input_s % 100) / 10;
    s3 = input_s % 10;
    cout << s1 + s2 + s3 << endl;
}
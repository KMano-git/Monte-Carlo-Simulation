#include <bits/stdc++.h>
using namespace std;

int pow(int x, int y) {
    int result = 1;
    for (int i = 0; i < y; i++) {
        result *= x;
    }
    return result;
}

int main() {
    int s;
    cin >> s;
    int count = 0;
    for (int i= 0; i < 3; i++) {
        if (s / pow(10, i) % 10 == 1) count++;
    }
    cout << count << endl;
}
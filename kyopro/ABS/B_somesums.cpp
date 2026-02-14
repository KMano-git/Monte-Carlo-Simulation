#include <bits/stdc++.h>
using namespace std;
int pow(int a, int b) {
    int res = 1;
    for (int i = 0; i < b; i++) {
        res *= a;
    }
    return res;
}

int main() {
    int n;
    cin >> n;
    int a, b;
    cin >> a >> b;
    
    int cnt = 0;
    for (int i = 1; i <= n; i++) {
        int sum = 0;
        for (int j = 0; j < 5; j++) sum += i / pow(10, j) % 10;
        if ( sum >= a && sum <= b ) cnt += i;
    }
    cout << cnt << endl;
}
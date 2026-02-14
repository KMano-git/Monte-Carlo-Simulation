#include <bits/stdc++.h>
using namespace std;

int main() {
    int a,b;
    cin >> a >> b;
    int count = 0;
    int n = b - a + 1; 
    for (int i = 0; i < n; i++) {
        int num = a + i;
        int s1, s2, s4, s5;
        s1 = num / 10000;
        s2 = (num / 1000) % 10;
        s4 = (num / 10) % 10;
        s5 = num % 10;
        if (s1 == s5 && s2 == s4) {
            count += 1;
        }
    }
    cout << count << endl;
}
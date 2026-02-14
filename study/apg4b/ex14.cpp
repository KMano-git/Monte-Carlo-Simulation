#include <bits/stdc++.h>
using namespace std;

int main() {
    int a, b, c;
    cin >> a >> b >> c;
    int heighest, lowest;

    heighest = max(a, max(b, c));
    lowest = min(a, min(b, c));

    cout << heighest - lowest << endl;
}
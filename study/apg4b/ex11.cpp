#include <bits/stdc++.h>
using namespace std;

int main() {
    int n, a;
    cin >> n >> a;
    int x = a;

    for (int i = 0; i < n; i++) {
        string op;
        int b;
        cin >> op >> b;
        if (op == "+") {
            x += b;
            cout << i + 1 << ":" << x << endl;
        } else if (op == "-") {
            x -= b;
            cout << i + 1 << ":" << x << endl;
        } else if (op == "*") {
            x *= b;
            cout << i + 1 << ":" << x << endl;
        } else if (op == "/" && b != 0) {
            // 0除算は許可しない
            x /= b;
            cout << i + 1 << ":" << x << endl;
        } else {
            cout << "error" << endl;
            break;
        }
    }
}
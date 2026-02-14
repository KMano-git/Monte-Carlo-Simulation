#include <bits/stdc++.h>
using namespace std;

int main() {
    int a, b;
    string op;
    cin >> a >> op >> b;
    if (op == "+") {
        cout << a + b << endl;
    } else if (op == "-") {
        cout << a - b << endl;
    } else if (op == "*") {
        cout << a * b << endl;
    } else if (op == "/" && b != 0) {
        // 0除算は許可しない
        cout << a / b << endl;
    } else {
        cout << "error" << endl;
    }
}
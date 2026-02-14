#include <bits/stdc++.h>
using namespace std;

int main() {
    int pattern;
    cin >> pattern;

    if (pattern == 1) {
        int price, quantity;
        cin >> price >> quantity;
        cout << price * quantity << endl;
    }
    else if (pattern == 2) {
        string text;
        int price, quantity;
        cin >> text >> price >> quantity;
        cout << text << "!" << endl;
        cout << price * quantity << endl;
    }

}
#include <bits/stdc++.h>
using namespace std;

int main() {
    int n, s;
    cin >> n >> s;
    vector<int> apple(n), pineapple(n);
    for (int i = 0; i < n; i++) {
        cin >> apple.at(i);
    }
    for (int i = 0; i < n; i++) {
        cin >> pineapple.at(i);
    }

    int count = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (apple.at(i) + pineapple.at(j) == s) {
                count++;
            }
        }
    }
    cout << count << endl;
}
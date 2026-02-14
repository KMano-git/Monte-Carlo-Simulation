#include <bits/stdc++.h>
using namespace std;

int main() {
    string S;
    cin >> S;

    int n = (S.size() - 1) / 2;
    int answer = 1;
    for (int i = 0; i < n; i++) {
        if (S.at(2 * i + 1) == '+') answer++;
        else answer--;
    }
    cout << answer << endl;
}
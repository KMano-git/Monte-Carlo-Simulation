#include <bits/stdc++.h>
using namespace std;

int main() {
    int N;
    cin >> N;
    set<string> user;
    for (int i = 0; i < N; i++) {
        string s;
        cin >> s;
        if (!user.count(s)){
            cout << i + 1 << endl;
            user.insert(s);
        }
    }
}
#include <bits/stdc++.h>
using namespace std;

int main() {
    string s;
    cin >> s;

    reverse(s.begin(), s.end());
    vector<string> words = { "maerd", "remaerd", "esare", "resare" };

    int pos = 0;
    string t = "";
    while (pos < s.size()) {
        bool found = false;
        for (int i = 0; i < words.size(); i++) {
            if (s.substr(pos, words[i].size()) == words[i]) {
                pos += words[i].size();
                t += words[i];                
                found = true;
            }
        }
        if (!found) {
            cout << "NO" << endl;
            return 0;
        }
    }
    if (t == s) cout << "YES" << endl;
    else cout << "NO" << endl;    
}
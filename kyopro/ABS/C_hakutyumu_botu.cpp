#include <bits/stdc++.h>
using namespace std;

int main() {
    string s;
    cin >> s;

    string t = ""; // dream, dreamer, erase, eraser
    int n = s.size();
    int cnt_d = 0;
    int cnt_s = 0;
    for (int i = 0; i < n; i++) {
        if (s.at(i) == 'd') cnt_d++;
        if (s.at(i) == 's') cnt_s++;
    }
    int cnt_er = 0;
    for (int i = 0; i < n - 1; i++) {
        if (s.at(i) == 'e' && s.at(i + 1) == 'r') cnt_er++, i++;
    }

    int cnt_dreamer = n        - 5*cnt_d - 4*cnt_s - cnt_er;
    int cnt_eraser  = 2*cnt_er + 5*cnt_d + 3*cnt_s - n;
    if (cnt_dreamer < 0 || cnt_dreamer > cnt_d) {
        cout << "NO" << endl;
        return 0;
    }else if (cnt_eraser < 0 || cnt_eraser > cnt_s) {
        cout << "NO" << endl;
        return 0;
    }
    
    // これだと、この検算方法を突破する特殊な状況が考えられるので、総当りで検証するコードが必要？
    // 総当りの書き方わからん
    // cnt_dream, cnt_dreamer, cnt_erase, cnt_eraserの値によっては4^cnt_word通りの計算が必要
    int cnt_dream = cnt_d - cnt_dreamer;
    int cnt_erase = cnt_s - cnt_eraser;
    int cnt_word = cnt_dream + cnt_dreamer + cnt_erase + cnt_eraser;
    vector<string> words = { "dream", "dreamer", "erase", "eraser" };
    
}
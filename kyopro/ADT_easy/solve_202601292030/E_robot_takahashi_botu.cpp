#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin >> n;
    string s;
    cin >> s;
    vector<int> W(n,0);
    int max_w = 0;
    for (int i = 0; i < n; i++) {
        cin >> W.at(i);
        max_w = max(max_w, W.at(i));
    }

    // 文字列sを数列a(n)に変換
    vector<int> a(n);
    for (int i = 0; i < n; i++) {
        a.at(i) = s.at(i) - '0';
    }

    // W(i)の値で区分けしたものが数列a(n)とどれだけ一致するかを出力する
    // 繰り返し処理において、インデックスで効率よく処理したい
    vector<int> W_sorted = W;
    // 0とmax_w+1を追加
    W_sorted.push_back(0);
    W_sorted.push_back(max_w + 1);
    sort(W_sorted.begin(), W_sorted.end());

    // 重複を削除
    W_sorted.erase(unique(W_sorted.begin(), W_sorted.end()), W_sorted.end());

    int max_count = 0;
    int W_sorted_size = W_sorted.size();
    for (int x = 0; x < W_sorted_size; x++) {
        int count = 0;
        vector<int> b(n); // x以上の値を1, 未満の値を0とする
        for (int i = 0; i < n; i++) {
            if (W.at(i) >= W_sorted.at(x)) {
                b.at(i) = 1;
            } else {
                b.at(i) = 0;
            }
        }
        for (int i = 0; i < n; i++) {
            if (a.at(i) == b.at(i)) {
                count++;
            }
        }
        max_count = max(max_count, count);
    }
    cout << max_count << endl;
}